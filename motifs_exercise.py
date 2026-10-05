from Bio import motifs
from Bio.Seq import Seq
import random
from Bio import SeqIO
import pandas as pd

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


class MotifFinder:
    def __init__(self, sequences, l, seed=None):
        """
        sequences: a list of strings
        l: length of the motif to find
        windows: a list with one list of l-mers per sequence: every sequence cut into all its l-mers once
        rnq: a random number generator with a seed for reproducibility
        """
        self.sequences = sequences
        self.l = l
        self.rng = random.Random(seed)  # a random number generator with a seed for reproducibility
        self.windows = [[sequence[i:i+self.l] for i in range(len(sequence) - self.l + 1)] for sequence in self.sequences]

    def total_distance(self, pattern):
        """
        pattern: a string of length l
        returns: int, like total_distance from Task 1, but using self.windows
        """
        total_distance = 0
        for window_list in self.windows:
            min_distance = self.l
            for window in window_list:
                distance = hamming_distance(pattern, window)
                if distance < min_distance:
                    min_distance = distance
            total_distance += min_distance
        return total_distance

    def median_string(self):
        """
        returns: the pattern with the smallest total distance out of all 4^l patterns of length l, together with that distance
        """
        best_pattern = None
        best_distance = self.l * len(self.sequences)  # maximum possible distance
        for i in range(4 ** self.l):
            pattern = ""
            for _ in range(self.l):
                pattern = "ACGT"[i % 4] + pattern
                i //= 4
            distance = self.total_distance(pattern)
            if distance < best_distance:
                best_distance = distance
                best_pattern = pattern
        return best_pattern, best_distance

    def randomized_search(self):
        """
        returns: the result of one run of randomized motif search
        """
        # pick a random l-mer in every sequence (self.rng.choice over its windows) → this is the current motif matrix;
        lmers = []
        for window in self.windows:
            lmers.append(self.rng.choice(window))
        while True:
            # build a MotifProfile from the current motifs (pseudocount 1)
            profile = MotifProfile(lmers, pseudocount=1)
            # in every sequence, take the profile-most-probable l-mer (most_probable_lmer) → new motifs
            new_motifs = []
            for sequence in self.sequences:
                new_motifs.append(profile.most_probable_lmer(sequence))

            # if the new motifs have a higher score (Task 1), keep them and go back to step 2; otherwise stop and return the current motifs and their score
            if score(new_motifs) > score(lmers):
                lmers = new_motifs
            else:
                return lmers, score(lmers)

    def best_of(self, runs):
        """
        runs: int, number of runs of randomized motif search
        returns: the best motifs and their score out of all runs
        """
        best_motifs = None
        best_score = 0
        for _ in range(runs):
            motifs, motifs_score = self.randomized_search()
            if motifs_score > best_score:
                best_score = motifs_score
                best_motifs = motifs
        return best_motifs, best_score

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


    # # check:
    # profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
    # print(profile.consensus())                       # ATGCGTA
    # print(round(profile.lmer_probability("ATGCGTA"), 4))  # 0.0122

    # two = MotifProfile(["GTAC", "TTAA"])
    # print(two.most_probable_lmer("ACTGGATGACCC"))    # TGAC
    # print(round(two.lmer_probability("TGAC"), 4))         # 0.0093

    # bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
    # bio.pseudocounts = 1
    # print(bio.consensus)     # ATGCGTA
    # print(bio.pwm["A"])      # the same numbers as your profile.ppm["A"]

    # rng = random.Random(1)                 # a random number generator with seed 1
    # i = rng.randint(0, 50)                 # random integer, 0 <= i <= 50 (both ends included!)
    # lmer = rng.choice(["ACG", "CGT", "GTA"])   # one random item of a list
    # print(i, lmer)

    # check on the lecture data:
    # lect_motif = MotifFinder(lecture_dna, 6)
    # print(lect_motif.median_string())  # "AGATAG",2
    # print(lect_motif.best_of(100))      # ['AGATAG', 'AGATAG', 'AGATAG', 'AGACAG', 'AGATAG', 'AGGTAG'], 34
    
    # find the planted motif:
    sequences = [str(record.seq) for record in SeqIO.parse("planted_motif.fasta", "fasta")]
    motif = MotifFinder(sequences, 7)

         # (['GCTACAG', 'GCTCAAG', 'ACTAAAG', 'GCTAATG', 'GATAAAG', 'GCTATAG', 'GATAAAG', 'TCTAAAG', 'GCTACAG', 'GCTAAGG'], 60)
    
    _, best_score = motif.best_of(200)
    # How many restarts do you need?
    median_pattern, _ = motif.median_string()
    N_trials = 50
    Rs = [5,10,20,50]
    successes = {R:{"measured_success_rate":0.0, "predicted_success_rate":0.0} for R in Rs}
    p_runs = 200
    ps = 0
    for _ in range(N_trials):
        if consensus(motif.randomized_search()[0]) == median_pattern:
            ps += 1
    p1 = ps / N_trials

    for R in Rs:
        succeed_count = 0
        for _ in range(N_trials):
            # print(f"best_of({R}): {motif.best_of(R)[1]}")
            if motif.best_of(R)[1] == best_score:
                succeed_count += 1
        successes[R]["measured_success_rate"] = succeed_count / N_trials
        successes[R]["predicted_success_rate"] = 1 - (1 - p1) ** R
    table = pd.DataFrame(successes)
    print(table)
    #                             5         10        20        50
    # measured_success_rate   0.220000  0.280000  0.500000  0.900000
    # predicted_success_rate  0.167498  0.306941  0.519669  0.840099

    # From the found motif to a scanner.
    best_motifs, _ = motif.best_of(100)
    found = motifs.create([Seq(m) for m in best_motifs])  # the motifs from best_of(100)
    found.pseudocounts = 1
    pssm = found.pssm

    n_windows = sum(len(w) for w in motif.windows)  # 10 * 54 = 540 per strand

    for fpr in (0.01, 0.001):
        threshold = pssm.distribution().threshold_fpr(fpr)
        total_hits = 0
        n_found = 0
        n_reverse = 0
        print(f"--- threshold_fpr({fpr}), score {threshold:.2f} ---")
        for idx, sequence in enumerate(sequences):
            for position, hit_score in pssm.search(Seq(sequence), threshold=threshold):
                total_hits += 1
                if position < 0:
                    # reverse strand: the window is on the forward strand at this place,
                    # and the motif matches its reverse complement
                    n_reverse += 1
                    start = len(sequence) + position
                else:
                    start = position
                    if sequence[start:start + motif.l] == best_motifs[idx]:
                        n_found += 1
                window = sequence[start:start + motif.l]
                strand = "-" if position < 0 else "+"
                print(f"seq {idx}: {strand} pos {start}, score {hit_score:.2f}, {window}")
        expected = 2 * n_windows * fpr  # both strands
        print(f"hits={total_hits}, found sites={n_found}, "
              f"reverse strand={n_reverse}, expected by chance={expected:.1f}")

        # expected results table:
        