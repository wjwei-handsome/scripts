
import sys

input_list_file = sys.argv[1]

with open(input_list_file) as f:
    data = f.read().splitlines()

# data = [1,2,3,20,25,28,104,114,124,5000,5005,8000]
data = [int(i) for i in data]
MAX_DIS = 200

res = []

for i in range(len(data)):
    if i == 0:
        res.append([data[i]])
    else:
        if data[i] - data[i-1] <= MAX_DIS:
            res[-1].append(data[i])
        else:
            res.append([data[i]])

for i in res:
    if len(i) > 1:
        print(",".join([str(j) for j in i])+'\n')
    else:
        continue

def countKmers(sequence: str, k: int) -> dict:
    """
    Gets the frequency of each k-mer in a string, skipping over non-base characters
    """
    KmerCount = {}

    baseCount = 0
    kmer = 0
    for i in range(len(sequence)):
        c = sequence[i]

        base = -1
        if c == 'A' or c == 'a':
            base = 0
        elif c == 'C' or c == 'c':
            base = 1
        elif c == 'G' or c == 'g':
            base = 2
        elif c == 'T' or c == 't':
            base = 3

        if base != -1:
            allButTwoHighest = ((1 << (2 * (k - 1))) - 1) & kmer
            kmer = (allButTwoHighest << 2) | base
            baseCount += 1

            if baseCount >= k:
                if kmer in KmerCount:
                    KmerCount[kmer] += 1
                else:
                    KmerCount[kmer] = 1
    return KmerCount

def jaccardSimilarity(seq1, seq2, k):
    sKmerFreq = countKmers(seq1, k)
    tKmerFreq = countKmers(seq2, k)

    intersection = 0
    union = 0

    for sKmer,sFreq in sKmerFreq.items():
        tFreq = tKmerFreq.get(sKmer, 0)
        intersection += min(sFreq, tFreq)
        union += max(sFreq, tFreq)

    for tKmer,tFreq in tKmerFreq.items():
        if tKmer not in sKmerFreq:
            union += tFreq

    return 1.0 * intersection / union

def editDistanceSimilarity(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    dp = [[0 for x in range(n + 1)] for x in range(m + 1)]
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0:
                dp[i][j] = j    # Min. operations = j
            elif j == 0:
                dp[i][j] = i    # Min. operations = i
            elif seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(dp[i][j-1],        # Insert
                                   dp[i-1][j],        # Remove
                                   dp[i-1][j-1])    # Replace
    return 1.0 * (m + n - dp[m][n]) / (m + n), 1.0 *(dp[n][m]/max(m,n))

def HammingDistanceSimilarity(seq1, seq2):
    distance = len([1 for c1, c2 in zip(seq1, seq2) if c1 == c2])
    return 1.0 * distance / len(seq1)