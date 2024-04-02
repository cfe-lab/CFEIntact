
from Bio import Align
from fractions import Fraction

aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -1.5
aligner.extend_gap_score = -0.2


class Alignment:
    def __init__(self, reference, query, score):
        self.reference = reference
        self.query = query
        self.score = score

    def __getitem__(self, index):
        if index == 0:
            return self.reference
        elif index == 1:
            return self.query
        else:
            raise IndexError("Index out of range")

    def distance(self):
        if self.score == 0:
            absolute = Fraction(0)
        else:
            absolute = Fraction(aligner.match_score - (self.score / len(self.reference)))

        # normalise it to the [0, 1) interval, with f(average_distance) = 0.5.
        average_distance: Fraction = Fraction("0.62")
        norm = Fraction(1) - (average_distance / (average_distance + absolute))

        return norm


def align(seq1, seq2):
    if seq1 and seq2:
        x = aligner.align(seq1, seq2)[0]
        return Alignment(x[0], x[1], x.score)
    elif seq1 or seq2:
        return align(seq1 or "-" * len(seq2), seq2 or "-" * len(seq1))
    else:
        return Alignment(seq1, seq2, 0)
