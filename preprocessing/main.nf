params.sraFile = "sra.tsv"

Channel.fromPath(params.sraFile)
        .splitCsv(sep: '\t') // [SRR924485, male_adipose]
        .set { sra_samples }

process test {
    tag "${prefix}"
    input:
    set val(sraID), val(sampleID) from sra_samples

    script:
    prefix = "${sraID}.${sampleID}"
    """
    echo "foo"
    """
}
