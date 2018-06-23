params.sraFile = "sra.tsv"
params.refFasta = "${params.ref_dir}/Gallus_gallus-5.0_genomic.fasta"
params.refAnno = "${params.ref_dir}/Gallus_gallus-5.0_genomic_cleaned.gff"

Channel.fromPath(params.refFasta).set { refFasta }
Channel.fromPath(params.refAnno).set { refAnno }
Channel.fromPath(params.sraFile)
        .splitCsv(sep: '\t') // [SRR924485, male, adipose]
        .set { sra_samples }

process index_ref {
    storeDir "${params.ref_dir}"

    input:
    file(fasta) from refFasta

    output:
    set file(fasta), file("${fasta}.*") into (indexed_ref, indexed_ref2)

    script:
    """
    echo "threads: \${NSLOTS:-\${threadsXL:-\${NTHREADS:-1}}}"

    hisat2-build "${fasta}" "${fasta}" -p \${NSLOTS:-\${threadsXL:-\${NTHREADS:-1}}}
    """
}

process align {
    tag "${prefix}"

    input:
    set val(sraID), val(sex), val(tissue), file(ref_fasta), file(ref_fasta_index: '*') from sra_samples.combine(indexed_ref)

    output:
    set val(sraID), val(sex), val(tissue), file("${output_bam}") into aligned_bams

    script:
    prefix = "${sraID}.${tissue}.${sex}"
    output_bam = "${prefix}.bam"
    """
    echo "threads: \${NSLOTS:-\${threadsMid:-\${NTHREADS:-1}}}"

    hisat2 -x ${ref_fasta} \
    --sra-acc ${sraID} \
    -p \${NSLOTS:-\${threadsMid:-\${NTHREADS:-1}}} | \
    samtools view -Sbo ${output_bam} -
    """
    // --rna-strandness $lib_type -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33
}

process sort {
    tag "${prefix}"
    publishDir "${params.output_dir}/alignments"

    input:
    set val(sraID), val(sex), val(tissue), file(bam) from aligned_bams

    output:
    set val(sraID), val(sex), val(tissue), file("${output_bam}") into sorted_bams

    script:
    prefix = "${sraID}.${tissue}.${sex}"
    output_bam = "${prefix}.sorted.bam"
    """
    echo "threads: \${NSLOTS:-\${threadsMid:-\${NTHREADS:-1}}}"

    sambamba sort \
    -t \${NSLOTS:-\${threadsMid:-\${NTHREADS:-1}}} \
    -o "${output_bam}" \
    ${bam}
    """
}

process assembly {
    tag "${prefix}"
    publishDir "${params.output_dir}"

    input:
    set val(sraID), val(sex), val(tissue), file(bam), file(anno) from sorted_bams.combine(refAnno)

    output:
    set val(sraID), val(sex), val(tissue), file("${output_gtf}") into assembled

    script:
    prefix = "${sraID}.${tissue}.${sex}"
    output_gtf = "${prefix}.gtf"
    """
    echo "threads: \${NSLOTS:-\${threadsMid:-\${NTHREADS:-1}}}"

    stringtie -G "${anno}" "${bam}" -o "${output_gtf}" -p \${NSLOTS:-\${threadsMid:-\${NTHREADS:-1}}}
    """
}

process filter {
    tag "${prefix}"
    publishDir "${params.output_dir}"

    input:
    set val(sraID), val(sex), val(tissue), file(gtf) from assembled

    output:
    set val(sraID), val(sex), val(tissue), file("${output_gtf}") into filtered_gtfs

    script:
    prefix = "${sraID}.${tissue}.${sex}"
    output_gtf = "${prefix}.filtered.gtf"
    """
    # make list of entries to remove from GTF:
    # remove 'transcript' without cov >2
    grep -w 'transcript' "${gtf}" | \
    grep -e 'cov "1\\.' -e 'cov "0\\.' | \
    cut -f9  | \
    cut -d " " -f 4 | \
    sort -u > remove.txt

    grep -vFf remove.txt "${gtf}" > "${output_gtf}"
    """
}

process extract_fasta {
    tag "${prefix}"
    publishDir "${params.output_dir}"

    input:
    set val(sraID), val(sex), val(tissue), file(gtf), file(ref_fasta), file(ref_fasta_index: '*') from filtered_gtfs.combine(indexed_ref2)

    output:
    set val(sraID), val(sex), val(tissue), file("${output_fasta}") into extracted_fastas

    script:
    prefix = "${sraID}.${tissue}.${sex}"
    output_fasta = "${prefix}.fasta"
    """
    gffread -w "${output_fasta}" -g "${ref_fasta}"
    """
}
// cut -f 9 $file | tr ";" "\n" | grep transcript | sed 's/^ //g' | cut -d " " -f 2 | grep -Ff - "$file" | gffread -w "${new}".fa -g "$referencegenome" - && python get_gene_length_filter.py "${new}".fa "${new}".putative_genes.fa && sed -i 's/ .*//' "${new}".putative_genes.fa
// get_gene_length_filter.py
#!/usr/bin/env python
//
// import sys
//
// infile = sys.argv[1] # input file
// outfile = sys.argv[2] # output file
//
// genes = {} # empty list
//
// with open(infile, "rU") as fh_in:
// 	for line in fh_in:
// 		line = line.strip()
// 		if line[0] == ">":
// 			gene_names = line
// 			genes[gene_names] = ''
// 		else:
// 			genes[gene_names]+=line
//
// for (name,val) in genes.items():
//     val = len(val)
//     with open(outfile, "a") as fh_out:
//         if val > 200:
//             fh_out.write(name)
//             fh_out.write("\n")
//             fh_out.write(genes[name])
//             fh_out.write("\n")


// bedtools intersect -wa -c -a SRR924202.bed -b chickspress_peptides.bed > srr924202_new_screen_intersect.txt
