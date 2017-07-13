#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <utility>
#include <algorithm>
#include <string>
#include <random>
#include <unordered_set>

using namespace std;

static size_t signatureWidth; // Signature size (in bits)
static size_t signatureSize;  // Signature size (in uint64_t)
static size_t kmerLength;     // Kmer length
static float density;         // % of sequence set as bits
static size_t clusterCount;   // # of clusters
static int kMeansIterations;  // # of iterations of k means
static bool fastaOutput;      // Output fasta or csv

vector<pair<string, string>> loadFasta(const char *path)
{
  vector<pair<string, string>> sequences;

  FILE *fp = fopen(path, "r");
  if (!fp) {
    fprintf(stderr, "Failed to load %s\n", path);
    exit(1);
  }
  for (;;) {
    char seqNameBuf[8192];
    if (fscanf(fp, " >%[^\n]\n", seqNameBuf) < 1) break;
    string sequenceBuf;

    for (;;) {
      int c = fgetc(fp);
      if (c == EOF || c == '>') {
        ungetc(c, fp);
        break;
      }
      if (isalpha(c)) {
        sequenceBuf.push_back(c);
      }
    }
    sequences.push_back(make_pair(string(seqNameBuf), sequenceBuf));
  }
  fclose(fp);

  return sequences;
}

void generateSignature(uint64_t *output, const pair<string, string> &fasta)
{
  // Generate a signature from the kmers contained within
  
  string fastaSequence = fasta.second;
  // If the sequence is shorter than the kmer length, pad it with Xs
  while (fastaSequence.size() < kmerLength) {
    fastaSequence.push_back('X');
  }
  
  ranlux24_base rng;
  uniform_int_distribution<int> dist(-64 * signatureSize, signatureSize * 64 - 1);
  vector<int> unflattenedSignature(signatureSize * 64);
  int setBits = density * signatureSize * 64;
  //fprintf(stderr, "%s\n", fastaSequence.c_str());
  
  for (size_t i = 0; i < fastaSequence.size() - kmerLength + 1; i++) {
    seed_seq rngSeed(begin(fastaSequence) + i, begin(fastaSequence) + i + kmerLength);
    rng.seed(rngSeed);
    string kmer(begin(fastaSequence) + i, begin(fastaSequence) + i + kmerLength);
    //fprintf(stderr, "- %s\n", kmer.c_str());
    
    for (int j = 0; j < setBits; j++) {
      int bitPos = dist(rng);
      if (bitPos >= 0) {
        unflattenedSignature[bitPos] += 1;
      } else {
        unflattenedSignature[bitPos + 64 * signatureSize] -= 1;
      }
    }
  }
  fill(output, output + signatureSize, 0);
  for (size_t i = 0; i < signatureSize * 64; i++) {
    if (unflattenedSignature[i] > 0) {
      output[i / 64] |= (uint64_t)1 << (i % 64);
    }
  }
}

vector<uint64_t> convertFastaToSignatures(const vector<pair<string, string>> &fasta)
{
  vector<uint64_t> output;
  // Allocate space for the strings
  
  output.resize(fasta.size() * signatureSize);
  
  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < fasta.size(); i++) {
    generateSignature(&output[signatureSize * i], fasta[i]);
  }
  
  return output;
}

vector<vector<size_t>> createClusterLists(const vector<size_t> &clusters)
{
  vector<vector<size_t>> clusterLists(clusterCount);
  for (size_t i = 0; i < clusters.size(); i++) {
    clusterLists[clusters[i]].push_back(i);
  }
  return clusterLists;
}

vector<uint64_t> createClusterSigs(const vector<vector<size_t>> &clusterLists, const vector<uint64_t> &sigs)
{
  vector<uint64_t> clusterSigs(signatureSize * clusterCount);
  #pragma omp parallel
  {
    vector<int> unflattenedSignature(signatureWidth);
    #pragma omp for
    for (size_t cluster = 0; cluster < clusterLists.size(); cluster++) {
      fill(begin(unflattenedSignature), end(unflattenedSignature), 0);
      
      for (size_t signature : clusterLists[cluster]) {
        const uint64_t *signatureData = &sigs[signatureSize * signature];
        for (size_t i = 0; i < signatureWidth; i++) {
          uint64_t signatureMask = (uint64_t)1 << (i % 64);
          if (signatureMask & signatureData[i / 64]) {
            unflattenedSignature[i] += 1;
          } else {
            unflattenedSignature[i] -= 1;
          }
        }
      }
      
      uint64_t *flattenedSignature = &clusterSigs[cluster * signatureSize];
      for (size_t i = 0; i < signatureWidth; i++) {
        if (unflattenedSignature[i] > 0) {
          flattenedSignature[i / 64] |= (uint64_t)1 << (i % 64);
        }
      }
    }
  }
  return clusterSigs;
}

void reclusterSignatures(vector<size_t> &clusters, const vector<uint64_t> &meanSigs, const vector<uint64_t> &sigs)
{
  #pragma omp parallel for
  for (size_t sig = 0; sig < clusters.size(); sig++) {
    const uint64_t *sourceSignature = &sigs[sig * signatureSize];
    size_t minHdCluster = 0;
    size_t minHd = numeric_limits<size_t>::max();
    for (size_t cluster = 0; cluster < clusterCount; cluster++) {
      const uint64_t *clusterSignature = &meanSigs[cluster * signatureSize];
      size_t hd = 0;
      for (size_t i = 0; i < signatureSize; i++) {
        hd += __builtin_popcountll(sourceSignature[i] ^ clusterSignature[i]);
      }
      if (hd < minHd) {
        minHd = hd;
        minHdCluster = cluster;
      }
    }
    clusters[sig] = minHdCluster;
  }
}

template<class RNG>
vector<uint64_t> createRandomSigs(RNG &&rng, const vector<uint64_t> &sigs)
{
  vector<uint64_t> clusterSigs(signatureSize * clusterCount);
  size_t signatureCount = sigs.size() / signatureSize;
  uniform_int_distribution<size_t> dist(0, signatureCount - 1);
  bool finished = false;
  
  unordered_set<string> uniqueSigs;
  for (size_t i = 0; i < signatureCount; i++) {
    size_t sig = dist(rng);
    string sigData(signatureSize * sizeof(uint64_t), ' ');
    memcpy(&sigData[0], &sigs[sig * signatureSize], signatureSize * sizeof(uint64_t));
    uniqueSigs.insert(sigData);
    if (uniqueSigs.size() >= clusterCount) {
      finished = true;
      break;
    }
  }
  
  if (!finished) {
    fprintf(stderr, "Error: there appear to be fewer unique signatures than clusters.\n");
    fprintf(stderr, "Reduce the cluster count\n");
    exit(1);
  }
  
  size_t i = 0;
  for (const auto &sig : uniqueSigs) {
    memcpy(&clusterSigs[i * signatureSize], sig.data(), signatureSize * sizeof(uint64_t));
    i++;
  }
  
  if (i != uniqueSigs.size()) {
    fprintf(stderr, "Mismatch: %zu != %zu\n", i, uniqueSigs.size());
    exit(1);
  }
  
  return clusterSigs;
}

vector<size_t> clusterSignatures(const vector<uint64_t> &sigs)
{
  auto rng = ranlux24_base();
  
  
  vector<size_t> clusters(sigs.size() / signatureSize);
  vector<vector<size_t>> clusterLists;
  vector<uint64_t> meanSigs;
  meanSigs = createRandomSigs(rng, sigs);
  for (int iteration = 0; iteration < kMeansIterations; iteration++) {
    fprintf(stderr, "Iteration %d\n", iteration);
    reclusterSignatures(clusters, meanSigs, sigs);
    clusterLists = createClusterLists(clusters);
    
    meanSigs = createClusterSigs(clusterLists, sigs);
  }
  
  size_t nonEmptyClusters = 0;
  for (const vector<size_t> &clusterList : clusterLists) {
    if (!clusterList.empty()) {
      nonEmptyClusters++;
    }
  }
  fprintf(stderr, "%llu/%llu non-empty clusters\n", static_cast<unsigned long long>(nonEmptyClusters), static_cast<unsigned long long>(clusterCount));
  
  return clusters;
}

void outputClusters(const vector<size_t> &clusters)
{
  for (size_t sig = 0; sig < clusters.size(); sig++)
  {
    printf("%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
  }
}

void outputFastaClusters(const vector<size_t> &clusters, const vector<pair<string, string>> &fasta)
{
  fprintf(stderr, "Writing out %zu records\n", clusters.size());
  for (size_t sig = 0; sig < clusters.size(); sig++)
  {
    printf(">%llu\n%s\n", static_cast<unsigned long long>(clusters[sig]), fasta[sig].second.c_str());
  }
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s (options) [fasta input]\n", argv[0]);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -sw [signature width]\n");
    fprintf(stderr, "  -k [kmer length]\n");
    fprintf(stderr, "  -d [signature density]\n");
    fprintf(stderr, "  -c [cluster count]\n");
    fprintf(stderr, "  -i [k-means iterations]\n");
    fprintf(stderr, "  --fasta-output\n");
    return 1;
  }
  signatureWidth = 256;
  kmerLength = 5;
  density = 1.0f / 21.0f;
  clusterCount = 1000;
  kMeansIterations = 4;
  fastaOutput = false;
  
  string fastaFile = "";
  
  for (int a = 1; a < argc; a++) {
    string arg(argv[a]);
    if (arg == "-sw") signatureWidth = atoi(argv[++a]);
    else if (arg == "-k") kmerLength = atoi(argv[++a]);
    else if (arg == "-d") density = atof(argv[++a]);
    else if (arg == "-c") clusterCount = atoi(argv[++a]);
    else if (arg == "-i") kMeansIterations = atoi(argv[++a]);
    else if (arg == "--fasta-output") fastaOutput = true;
    else if (fastaFile.empty()) fastaFile = arg;
    else {
      fprintf(stderr, "Invalid or extra argument: %s\n", arg.c_str());
      exit(1);
    }
  }
    
  if (signatureWidth <= 0 || signatureWidth % 64 != 0) {
    fprintf(stderr, "Error: signature width is not a multiple of 64\n");
    return 1;
  }
  if (kmerLength <= 0) {
    fprintf(stderr, "Error: kmer length must be a positive nonzero integer\n");
    return 1;
  }
  if (density < 0.0f || density > 1.0f) {
    fprintf(stderr, "Error: density must be a positive value between 0 and 1\n");
    return 1;
  }
  if (clusterCount <= 0) {
    fprintf(stderr, "Error: cluster count must be a positive nonzero integer\n");
    return 1;
  }
  if (kMeansIterations <= 0) {
    fprintf(stderr, "Error: iterations must be a positive nonzero integer\n");
    return 1;
  }
  
  signatureSize = signatureWidth / 64;
  
  fprintf(stderr, "Loading fasta...");
  auto fasta = loadFasta(fastaFile.c_str());
  fprintf(stderr, " loaded %llu sequences\n", static_cast<unsigned long long>(fasta.size()));
  fprintf(stderr, "Converting fasta to signatures...");
  auto sigs = convertFastaToSignatures(fasta);
  fprintf(stderr, " done\n");
  fprintf(stderr, "Clustering signatures...\n");
  auto clusters = clusterSignatures(sigs);
  fprintf(stderr, "Writing output\n");
  if (!fastaOutput) {
    outputClusters(clusters);
  } else {
    outputFastaClusters(clusters, fasta);
  }
  
  return 0;
}
