params.sraFile = "sra.tsv"
params.refFasta = "${params.ref_dir}/Gallus_gallus-5.0_genomic.fasta"

Channel.fromPath(params.refFasta).set { refFasta }

Channel.fromPath(params.sraFile)
        .splitCsv(sep: '\t') // [SRR924485, male_adipose]
        .set { sra_samples }

process index_ref {
    storeDir "${params.ref_dir}"

    input:
    file(fasta) from refFasta

    output:
    set file(fasta), file("${fasta}.*") into indexed_ref

    script:
    """
    echo "threads: \${NSLOTS:-\${threadsXL:-\${NTHREADS:-1}}}"
    hisat2-build "${fasta}" "${fasta}" -p \${NSLOTS:-\${threadsXL:-\${NTHREADS:-1}}}
    """
}
// indexed_ref.toList().set { indexed_ref_all }
// sra_samples.combine(indexed_ref).subscribe { println "[sra_samples] ${it}" }

process align {
    tag "${prefix}"

    input:
    set val(sraID), val(sampleID), file(ref_fasta), file(ref_fasta_index: '*') from sra_samples.combine(indexed_ref)

    script:
    prefix = "${sraID}.${sampleID}"
    output_bam = "${prefix}.bam"
    """
    echo "threads: \${NSLOTS:-\${threadsMid:-\${NTHREADS:-1}}}"

    hisat2 -x ${ref_fasta} \
    --sra-acc ${sraID} \
    -p \${NSLOTS:-\${threadsMid:-\${NTHREADS:-1}}} | \
    samtools view -Sbo ${output_bam} -
    """
        // --rna-strandness $lib_type \
    // -5 $five_trim \
    // -3 $three_trim \
    // --min-intronlen $min_intl \
    // --max-intronlen $max_intl \
    // --dta --phred33
}
