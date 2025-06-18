process EXTRACT_CGC_RANGES {
    input:
    path cgc_standard_out
    path substrate_prediction

    output:
    path "cgc_ranges.tsv"

    script:
    """
    python3 ../../bin/extract_cgc_ranges.py ${cgc_standard_out} ${substrate_prediction} cgc_ranges.tsv
    """
}
