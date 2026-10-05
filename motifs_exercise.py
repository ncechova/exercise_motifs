from Bio import motifs
from Bio.Seq import Seq

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
    """
    motifs: a list of equally long strings
    returns: a dictionary with counts of each nucleotide at each position
    """
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
    """
    motifs: a list of equally long strings
    returns: the score of the motifs, which is the sum of the counts of the most common nucleotide at each position
    """
    counts = count_matrix(motifs)
    score = 0
    for i in range(len(motifs[0])):
        max_count = max(counts[nucleotide][i] for nucleotide in 'ACGT')
        score += max_count
    return score

def consensus(motifs):
    """
    motifs: a list of equally long strings
    returns: the consensus sequence of the motifs
    """
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
    """
    a, b: two strings of equal length
    returns: the Hamming distance between a and b
    """
    distance = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            distance += 1
    return distance

def total_distance(pattern, sequences):
    """
    pattern: a string
    sequences: a list of strings
    returns: the total distance between the pattern and the sequences, which is the sum of the
    minimum Hamming distances between the pattern and each sequence
    """
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

# task 2
class MotifProfile:
    def __init__(self, motifs, pseudocount=1):
        """
        motifs: a list of equally long strings
        pseudocount: added to every count to avoid zero probabilities
        """
        self.motifs = motifs
        self.pseudocount = pseudocount
        self.l = len(motifs[0])
        self.t = len(motifs)
        self.ppm = self._calculate_ppm()

    def _calculate_ppm(self):
        counts = count_matrix(self.motifs)
        ppm = {
        'A': [0.0] * self.l,
        'C': [0.0] * self.l,
        'G': [0.0] * self.l,
        'T': [0.0] * self.l
        }
        for nucleotide in 'ACGT':
             for i, count in enumerate(counts[nucleotide]):
                ppm[nucleotide][i] += (count + self.pseudocount) / (self.t + 4 * self.pseudocount) 
        return ppm

    def lmer_probability(self, lmer):
        """
        lmer: a string of length l
        returns: float, the product of the probabilities of its bases, column by column
        """
        probability = 1.0
        for i, nucleotide in enumerate(lmer):
            probability *= self.ppm[nucleotide][i]
        return probability

    def most_probable_lmer(self, sequence):
        """
        sequence: a string at least self.l long
        returns: the lmer of length l in the sequence with the highest probability
        """
        max_probability = 0.0
        most_probable_lmer = ""
        for i in range(len(sequence) - self.l + 1):
            lmer = sequence[i:i+self.l]
            probability = self.lmer_probability(lmer)
            if probability > max_probability:
                max_probability = probability
                most_probable_lmer = lmer
        return most_probable_lmer

    def consensus(self):
        """
        returns: str of length self.l, the most probable base of every column
        """
        consensus = ""
        for i in range(self.l):
            max_prob = 0
            best_nucleotide = ''
            for nucleotide in 'ACGT':
                if self.ppm[nucleotide][i] > max_prob:
                    max_prob = self.ppm[nucleotide][i]
                    best_nucleotide = nucleotide
            consensus += best_nucleotide
        return consensus
 


if __name__ == "__main__":
    # counts = count_matrix(lecture_dna)
    # print(count_matrix(["ACGT", "ATGT", "CCGA"]))
    # print(score(["ACGT", "ATGT", "CCGA"]))
    # print(consensus(["AC", "GT"]))
    # print(hamming_distance("ACGT", "ACCA"))
    # print(total_distance("AC", ["GACT", "TTAG"]))


    # red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]  # slide 19
    # best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]

    # print(score(red))                                # 26
    # print(consensus(best), score(best))              # AGATAG 34
    # print(hamming_distance("TAAGTT", "TGAATT"))      # 2
    # print(total_distance("TGCGTT", lecture_dna))     # 13

    # profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
    # print(profile.l)           # 7
    # print(profile.ppm["A"])    # [0.5, 0.25, 0.125, 0.125, 0.25, 0.125, 0.5]
    # print(profile.lmer_probability("ATGCGTA")) # 0.0122
    # print(MotifProfile(["GTAC", "TTAA"]).most_probable_lmer("ACTGGATGACCC")) #"TGAC"
    # print(MotifProfile(["GTAC", "TTAA"]).most_probable_lmer("GGATGA")) # "GGAT"
    # print(profile.consensus()) # "ATGCGTA"


    # check:
    profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
    print(profile.consensus())                       # ATGCGTA
    print(round(profile.lmer_probability("ATGCGTA"), 4))  # 0.0122

    two = MotifProfile(["GTAC", "TTAA"])
    print(two.most_probable_lmer("ACTGGATGACCC"))    # TGAC
    print(round(two.lmer_probability("TGAC"), 4))         # 0.0093

    bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
    bio.pseudocounts = 1
    print(bio.consensus)     # ATGCGTA
    print(bio.pwm["A"])      # the same numbers as your profile.ppm["A"]
    print(profile.ppm["A"]) 
