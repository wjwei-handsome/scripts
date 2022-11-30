
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

from itertools import groupby

def cigarCategory(alignmentColumn):
    x, y = alignmentColumn
    if x == "-":
        if y == "-":
            return "P"
        else:
            return "I"
    else:
        if y == "-":
            return "D"
        elif x != y:
            return "X"
        else:
            return "M"

def cigarParts(beg, alignmentColumns, end):
    if beg:
        yield str(beg) + "H"
    # (doesn't handle translated alignments)
    for k, v in groupby(alignmentColumns, cigarCategory):
        yield str(sum(1 for _ in v)) + k
    if end:
        yield str(end) + "H"

def get_cigar(m1, m2):
    # qRevStart = m2.sequence_size - m2.alignment_start - m2.alignment_size
    cigar = "".join(cigarParts(0, zip(m1, m2), 0))
    return cigar

from itertools import combinations
from operator import mul
from functools import reduce
import random

def get_random_01_list(k: int) -> list:
    """
    Get a random list of 0 and 1
    """
    return [random.randint(0, 1) for _ in range(k)]

test_data = {i:get_random_01_list(5) for i in range(1,10)}

for i in combinations(test_data.keys(), 3):
    a=test_data[i[0]]
    b=test_data[i[1]]
    c=test_data[i[2]]
    print(f"{i}:{reduce(mul,a)}\t{reduce(mul,b)}\t{reduce(mul,c)}")
