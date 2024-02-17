# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    score = -8
    A = "MQ-R"
    B = "MYQR"
    score = -8
    nw = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    alignment_score = nw.align(seq1, seq2)

    assert alignment_score[0] == score,"Incorrected score"
    assert alignment_score[1] == B,"Incorrect alignment"
    assert alignment_score[2] == A,"Incorrect alignment"
    assert np.all(nw._align_matrix == np.array ([[  0, np.NINF, np.NINF, np.NINF,np.NINF],[np.NINF,1,-11,-12,-13],[np.NINF,-11,1,-9,-11],[np.NINF,-12,-10,1,-8]])),"Incorrect score matrix"
    assert np.all(nw._gapB_matrix == np.array([[[0,np.NINF,np.NINF,np.NINF,np.NINF],[-11,np.NINF,np.NINF,np.NINF,np.NINF],[-12, -10, -22, -23, -24],[-13, -11, -10, -20, -22]]])),"Incorrect Ix matrix"
    assert np.all(nw._gapA_matrix == np.array([[ 0, -11, -12, -13, -14,],[np.NINF,np.NINF, -10, -11, -12],[np.NINF,np.NINF, -22, -10, -11],[np.NINF,np.NINF, -23, -21, -10,]])),"Incorrect Iy matrix"

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    alignment_score = nw.align(seq3, seq4)
    assert np.all(nw._back == np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],[0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2],[0, 1, 0, 2, 2, 2, 2, 2, 2, 2, 2],[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]])),"Incorrect backtrace matrix"
    assert alignment_score[0] == -7,"Incorrect alignment score"
    assert alignment_score[1] == "MAVHQLIRRP","Incorrect alignment"
    assert alignment_score[2] == "MQ---LIRHP","Incorrect alignment"


