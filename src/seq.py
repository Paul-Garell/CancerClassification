# pip install biopython
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
import Ref.TP53Ref as TP53Ref
import Ref.GATA3Ref as GATA3Ref


def align(input):
    # Aligment Parameters
    matrix = substitution_matrices.load("BLOSUM62")
    gap_open = -10
    gap_extend = -0.5

    max_align, max_gene = None, None

    # TP53 Gene Aligment
    align_TP53 = pairwise2.align.localds(
        TP53Ref.TP53Ref, input, match_dict=matrix, open=gap_open, extend=gap_extend
    )[0]
    max_align, max_gene = align_TP53.score, "TP53"

    # GATA3 Gene Aligment
    align_GATA3 = pairwise2.align.localds(
        GATA3Ref.GATA3Ref, input, match_dict=matrix, open=gap_open, extend=gap_extend
    )[0]
    if align_GATA3.score > max_align:
        max_align, max_gene = align_GATA3.score, "GATA3"

    print("Score: " + str(max_align))
    print("Gene: " + max_gene)
    return max_gene

    # Utilize gene frequency (mult score by percentage) -> choose max ->
    # Gene (TP53), Start, End, Ref, Var
