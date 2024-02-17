# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    nw = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/PAM100.mat", gap_open=-10, gap_extend=-1)
    alignment = nw.align(hs_seq,gg_seq)
    alignment2 = nw.align(hs_seq,mm_seq)
    alignment3 = nw.align(hs_seq,br_seq)
    alignment4 = nw.align(hs_seq,tt_seq)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    list_names = ["Gallus_gallus","Mus_musculus","Balaeniceps_rex","tursiops_truncatus"] 
    scores = [alignment[0],alignment2[0],alignment3[0],alignment4[0]]
    dict_scores = {}
    for i in range(len(list_names)):
        dict_scores[list_names[i]] = scores[i]
    sorted_dict = dict(sorted(dict_scores.items(), key=lambda item: item[1],reverse=True))
    print("Original:",dict_scores)
    print("Sorted (descending order):",sorted_dict)
    
    

if __name__ == "__main__":
    main()
