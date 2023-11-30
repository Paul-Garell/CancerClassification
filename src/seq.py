# pip install biopython
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
import Ref.TP53Ref as TP53Ref
import Ref.GATA3Ref as GATA3Ref
import Ref.PIK3CARef as PIK3CARef


def get_max(input):
    # Aligment Parameters
    matrix = substitution_matrices.load("BLOSUM62")
    gap_open = -10
    gap_extend = -0.5

    # TODO: still need to do CDH1
    max_align, max_gene, max_res, ref = -float("inf"), None, None, None

    gene_ref = {
        "TP53": TP53Ref.TP53Ref,
        "GATA3": GATA3Ref.GATA3Ref,
        "PIK3CA": PIK3CARef.PIK3CARef,
    }

    for gene in gene_ref:
        align = pairwise2.align.localds(
            gene_ref[gene],
            input,
            match_dict=matrix,
            open=gap_open,
            extend=gap_extend,
        )[0]
        if align.score > max_align:
            max_align, max_gene, max_res = align.score, gene, align
            ref = gene_ref[gene]

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
