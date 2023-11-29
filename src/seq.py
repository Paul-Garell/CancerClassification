# pip install biopython
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
import Ref.TP53Ref as TP53Ref
import Ref.GATA3Ref as GATA3Ref
import Ref.PIK3CARef as PIK3CARef


def get_max(input):
    simple = True

    # Aligment Parameters
    matrix = substitution_matrices.load("BLOSUM62")
    gap_open = -10
    gap_extend = -0.5

    max_align, max_gene, max_res, ref = None, None, None, None

    # TP53 Gene Aligment
    align_TP53 = pairwise2.align.localds(
        TP53Ref.TP53Ref, input, match_dict=matrix, open=gap_open, extend=gap_extend
    )[0]
    max_align, max_gene, max_res = align_TP53.score, "TP53", align_TP53
    ref = TP53Ref.TP53Ref

    # GATA3 Gene Aligment
    align_GATA3 = pairwise2.align.localds(
        GATA3Ref.GATA3Ref, input, match_dict=matrix, open=gap_open, extend=gap_extend
    )[0]
    if align_GATA3.score > max_align:
        max_align, max_gene, max_res = align_GATA3.score, "GATA3", align_GATA3
        ref = GATA3Ref.GATA3Ref

    if simple:
        return (max_gene, max_res, ref)

    # PIK3CA Gene Aligment
    align_PIK3CA = pairwise2.align.localds(
        PIK3CARef.PIK3CARef, input, match_dict=matrix, open=gap_open, extend=gap_extend
    )[0]
    if align_PIK3CA.score > max_align:
        max_align, max_gene, max_res = align_PIK3CA.score, "PIK3CA", align_PIK3CA
        ref = PIK3CARef.PIK3CARef

    # TODO: still need to do CDH1

    return (max_gene, max_res, ref)


def align(input):
    # TODO: Eventually we will utilize gene freq (mult score by percentage)
    start_pos = {
        "TP53": 7571720,
        "GATA3": 8096667,
        "PIK3CA": 178866311,
        "CDH1": 68771195,
    }

    max_gene_name, max_res, ref = get_max(input)
    result_sequence = max_res.seqB[max_res.start : max_res.end]

    # Compares result sequence vs reference seqeuence
    difference = []
    for i in range(0, max_res.end - max_res.start):
        ref_base = ref[i + max_res.start]
        res_base = result_sequence[i]
        if ref_base != res_base:
            difference.append((ref_base, res_base))

    # [Gene, Start, End, [(Ref, Var)]]
    # Ex. ["TP53", 8106058, 8106058, [(T, A)]]
    res = [
        max_gene_name,
        max_res.start + start_pos[max_gene_name],
        max_res.end + start_pos[max_gene_name],
        difference,
    ]

    print("Res Sequence: " + result_sequence)
    print("Ref Sequence: " + ref[max_res.start : max_res.end])
    print("Result Vector: ")
    print(res)
    print()

    return res
