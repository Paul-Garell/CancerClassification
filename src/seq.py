# pip install biopython
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
import Ref.TP53Ref as TP53Ref
import Ref.GATA3Ref as GATA3Ref
import Ref.PIK3CARef as PIK3CARef
import Ref.CDH1Ref as CDH1Ref


def get_max(input):
    # Aligment Parameters
    matrix = substitution_matrices.load("BLOSUM62")
    gap_open = -10
    gap_extend = -0.5

    max_align, max_gene, max_res, ref = -float("inf"), None, None, None

    gene_ref = {
        "TP53": TP53Ref.TP53Ref,
        "GATA3": GATA3Ref.GATA3Ref,
        # "PIK3CA": PIK3CARef.PIK3CARef,
        # "CDH1": CDH1Ref.CDH1Ref,
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
    info = {
        "TP53": [7571720, 17],
        "GATA3": [8096667, 10],
        "PIK3CA": [178866311, 3],
        "CDH1": [68771195, 16],
    }

    max_gene_name, max_res, ref = get_max(input)
    result_sequence = max_res.seqB[max_res.start : max_res.end]

    # Compares result sequence vs reference seqeuence
    difference = []
    for i in range(0, max_res.end - max_res.start):
        ref_base = ref[i + max_res.start]
        res_base = result_sequence[i]
        if ref_base != res_base:
            difference.append(
                [ref_base, res_base, i + max_res.start + info[max_gene_name][0]]
            )

    res = []
    for dif in difference:
        res.append(
            [
                info[max_gene_name][1],  # chromosome
                dif[2],  # start
                dif[2],  # end
                dif[0],  # ref
                dif[1],  # variation
                max_gene_name,  # name
            ]
        )
    print(res)

    return res
