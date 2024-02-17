# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        #converting strings to list for ease: i added this, this was not here before
        seqA1 = [i for i in seqA]
        seqB1 = [i for i in seqB]

        seqA1.insert(0," ")
        seqB1.insert(0," ")
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing

        #score matrices
        Ix = np.zeros((len(seqB1),len(seqA1))) #score of alignment if seqA aligned w gap #i swapped A and B
        Iy = np.zeros((len(seqB1),len(seqA1))) #score of alignment if seqB aligned with gap
        M = np.zeros((len(seqB1),len(seqA1))) #score of alignment matrix

        #traceback matrices
        M_back = np.zeros((len(seqB1),len(seqA1)),dtype=object) #removed +1 from all the lengths because i forgot i accounted for that by adding a space as the first element of the sequence lists
        Iy_back = np.zeros((len(seqB1),len(seqA1)),dtype=object)
        Ix_back = np.zeros((len(seqB1),len(seqA1)),dtype=object)

        for i in range(1,len(seqB1)):
            #M[i+1][0] = self.gap_open + (i * self.gap_extend)
            M[i][0] = np.NINF #removed [i+1]
            Iy[i][0] = np.NINF
            Ix[i][0] = self.gap_open + (i * self.gap_extend) #this was also np.NINF

        for j in range(1,len(seqA1)):
            M[0][j] = np.NINF
            #M[0][j+1] = self.gap_open + (j * self.gap_extend)
            Ix[0][j] = np.NINF
            Iy[0][j] = self.gap_open + (j * self.gap_extend)
        #Ix[0][0] = self.gap_open
        #Iy[0][0] = self.gap_open
        M[0][0] = 0

        # TODO: Implement global alignment here
        def score(x,y):
            if x==y:
                return 1
            else:
                return 0 #video uses -1 penalty
        for i in range(1,len(seqB1)):
            for j in range(1,len(seqA1)):
                Iy[i][j] = max(M[i][j-1] + self.gap_open + self.gap_extend, Iy[i][j-1] + self.gap_extend)  #swapped x and y here
                Ix[i][j] = max(M[i-1][j] + self.gap_open + self.gap_extend, Ix[i-1][j] + self.gap_extend)
                M[i][j] = max((M[i-1][j-1] + score(seqA1[j], seqB1[i])), (Ix[i-1][j-1] + score(seqA1[j], seqB1[i])), (Iy[i-1][j-1] + score(seqA1[j], seqB1[i]))) #M[i-1][j]+self.gap_open, M[i][j-1]+self.gap_open) 

                #backtrace matrices
                Ix_back[i][j] = 0 if M[i-1][j] + self.gap_open + self.gap_extend >= Ix[i-1][j] + self.gap_extend else 1
                Iy_back[i][j] = 0 if M[i][j-1] + self.gap_open + self.gap_extend >= Iy[i][j-1] + self.gap_extend else 1

                if M[i][j] == M[i-1][j-1] + score(seqA1[j], seqB1[i]):
                    M_back[i][j] = 0 #diagonal
                elif M[i][j] == Ix[i-1][j-1] + score(seqA1[j], seqB1[i]):
                    M_back[i][j] = 1 #came from Ix
                elif M[i][j] == Iy[i-1][j-1] + score(seqA1[j], seqB1[i]):
                    M_back[i][j] = 2 #came from Iy
        
        self._align_matrix = M
        self._gapA_matrix = Iy
        self._gapB_matrix = Ix

        self._back = M_back
        self._back_A = Iy_back
        self._back_B = Ix_back

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        
        i, j = len(self._seqB), len(self._seqA)  
        aligned_seqA = []
        aligned_seqB = []

        #Determine the matrix from which to start: it could be M, Ix, or Iy
        current_matrix = 'M' if self._align_matrix[i][j] >= self._gapB_matrix[i][j] and self._align_matrix[i][j] >= self._gapA_matrix[i][j] else ('Ix' if self._gapB_matrix[i][j] > self._gapA_matrix[i][j] else 'Iy') #max number
        while i > 0 or j > 0:
            if current_matrix == 'M':
                if self._back[i][j] == 0:  #Diagonal
                    aligned_seqA.append(self._seqA[j-1])
                    aligned_seqB.append(self._seqB[i-1])
                    i -= 1
                    j -= 1
                elif self._back[i][j] == 1:  #Came from Ix
                    current_matrix = 'Ix'
                else:  #Came from Iy
                    current_matrix = 'Iy'
            elif current_matrix == 'Ix':
                if self._back_B[i][j] == 0:  #Up (gap in seqA) #used to be 0
                    aligned_seqA.append('-')
                    aligned_seqB.append(self._seqB[i-1])
                    i -= 1
                #Switch to M if this was a gap opening
                    if i > 0 and self._align_matrix[i][j] + self.gap_open + self.gap_extend >= self._gapB_matrix[i][j] + self.gap_extend:
                        current_matrix = 'M'
                else:  #Continue in Ix
                    aligned_seqA.append('-')
                    aligned_seqB.append(self._seqB[i-1])
                    i -= 1
            elif current_matrix == 'Iy':
                if self._back_A[i][j] == 0:  #Left (gap in seqB) #used to be 0
                    aligned_seqA.append(self._seqA[j-1])
                    aligned_seqB.append('-')
                    j -= 1
                #Switch to M if this was a gap opening
                    if j > 0 and self._align_matrix[i][j] + self.gap_open + self.gap_extend >= self._gapA_matrix[i][j] + self.gap_extend:
                        current_matrix = 'M'
                else:  #In Iy
                    aligned_seqA.append(self._seqA[j-1])
                    aligned_seqB.append('-')
                    j -= 1

    #Reverse the alignments 
        self.seqA_align = ''.join(reversed(aligned_seqA))
        self.seqB_align = ''.join(reversed(aligned_seqB))

    #Compute alignment score 
        i, j = len(self._seqB), len(self._seqA)
        self.alignment_score = max(self._align_matrix[i][j], self._gapB_matrix[i][j], self._gapA_matrix[i][j])  #Choose the maximum score from the three matrices

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header

#nw = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

# Assuming you have two sequences seqA and seqB to align
#seqA = "AAT"
#seqB = "ACACT"

# Call the align method with the two sequences
#alignment_score = nw.align(seqA, seqB)
#print(alignment_score)