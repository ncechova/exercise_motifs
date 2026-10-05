lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]

# task 1
def count_matrix(motifs):
    counts = {
        'A': [0] * len(motifs[0]),
        'C': [0] * len(motifs[0]),
        'G': [0] * len(motifs[0]),
        'T': [0] * len(motifs[0])
    }
    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            counts[nucleotide][i] += 1
    return counts

def score(motifs):
    counts = count_matrix(motifs)
    score = 0
    for i in range(len(motifs[0])):
        max_count = max(counts[nucleotide][i] for nucleotide in 'ACGT')
        score += max_count
    return score

def consensus(motifs):
    counts = count_matrix(motifs)
    consensus = ""
    for i in range(len(motifs[0])):
        max_count = 0
        best_nucleotide = ''
        for nucleotide in 'ACGT':
            if counts[nucleotide][i] > max_count:
                max_count = counts[nucleotide][i]
                best_nucleotide = nucleotide
        consensus += best_nucleotide
    return consensus


def hamming_distance(a, b):
    distance = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            distance += 1
    return distance

def total_distance(pattern, sequences):
    l = len(pattern)
    total_distance = 0
    for sequence in sequences:
        min_distance = l
        for i in range(len(sequence) - l + 1):
            window = sequence[i:i+l]
            distance = hamming_distance(pattern, window)
            if distance < min_distance:
                min_distance = distance
        total_distance += min_distance
    return total_distance


if __name__ == "__main__":
    # counts = count_matrix(lecture_dna)
    # print(count_matrix(["ACGT", "ATGT", "CCGA"]))
    # print(score(["ACGT", "ATGT", "CCGA"]))
    # print(consensus(["AC", "GT"]))
    # print(hamming_distance("ACGT", "ACCA"))
    # print(total_distance("AC", ["GACT", "TTAG"]))
    red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]  # slide 19
    best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]

    print(score(red))                                # 26
    print(consensus(best), score(best))              # AGATAG 34
    print(hamming_distance("TAAGTT", "TGAATT"))      # 2
    print(total_distance("TGCGTT", lecture_dna))     # 13