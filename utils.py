CHROMOSOME_ORDER = [
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY"
]

def sort_by_chrom(x):
    """
    Sort a list by the predefined chromosome order

    Parameters
    ----------
    x: tuple (chromosome, size)
    """
    return CHROMOSOME_ORDER.index(x[0]) if x[0] in CHROMOSOME_ORDER else 999
