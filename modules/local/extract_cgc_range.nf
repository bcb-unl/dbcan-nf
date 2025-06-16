process EXTRACT_CGC_RANGES {
    input:
    path cgc_tsv
    path substrate_tsv

    output:
    path "cgc_ranges.tsv"

    script:
    """
    python3 scripts/extract_cgc_ranges.py ${cgc_tsv} ${substrate_tsv} cgc_ranges.tsv
    """
}
