params.sraFile = "sra.tsv"
params.refFasta = "${params.ref_dir}/Gallus_gallus-5.0_genomic.fasta"
params.refAnno = "${params.ref_dir}/Gallus_gallus-5.0_genomic_cleaned.gff"
params.proteomeFastaGz = "${params.ref_dir}/uniprot-proteome_gallus_gallus_AUP000000539.fasta.gz"
// http://www.uniprot.org/uniprot/?query=gallus+gallus&sort=score
// http://www.uniprot.org/help/programmatic_access
params.blast_index_basename = "blast.uniprot-chicken"

Channel.fromPath(params.refFasta).set { refFasta }
Channel.fromPath(params.refAnno).set { refAnno }
Channel.fromPath(params.proteomeFastaGz).set { proteomeFastaGz }
// Channel.fromPath(params.sraFile)
//         .splitCsv(sep: '\t') // [SRR924485, male, adipose]
//         .set { sra_samples }
Channel.from([
    ["SRR924485", "male", "adipose"],
    ["SRR924539", "female", "adipose"]
    ])
    .set { sra_samples }

process unzip_proteomeFastaGz {
    storeDir "${params.ref_dir}"

    input:
    file(fastaGz) from proteomeFastaGz

    output:
    file("${output_fasta}") into proteomeFasta

    script:
    output_fasta = "${fastaGz}".replaceFirst(/.gz$/, "")
    """
    gunzip -c "${fastaGz}" > "${output_fasta}"
    """
}

process blast_index {
    storeDir "${params.ref_dir}"

    input:
    file(fasta) from proteomeFasta

    output:
    file("${params.blast_index_basename}*") into blast_index_files

    script:
    """
    makeblastdb -in "${fasta}" -dbtype prot -out "${params.blast_index_basename}"
    """
}

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
    gffread "${gtf}" -w "${output_fasta}" -g "${ref_fasta}"
    """
}

process filter_fasta {
    tag "${prefix}"
    publishDir "${params.output_dir}"

    input:
    set val(sraID), val(sex), val(tissue), file(fasta) from extracted_fastas

    output:
    set val(sraID), val(sex), val(tissue), file("${output_fasta}") into filtered_fastas

    script:
    prefix = "${sraID}.${tissue}.${sex}"
    output_fasta = "${prefix}.filtered.fasta"
    """
    #!/usr/bin/env python
    from Bio import SeqIO

    input_fasta = "${fasta}"
    output_fasta = "${output_fasta}"

    def filter_fasta(input_fasta, seq_len = 200):
        '''
        Emit only .fasta records with a sequence length greater than the limit
        '''
        with open(input_fasta, "rU") as fasta_in:
            for record in SeqIO.parse(fasta_in, "fasta"):
                if len(record.seq) > seq_len:
                    yield(record)

    SeqIO.write(filter_fasta(input_fasta = input_fasta), output_fasta, "fasta")
    """
}

blast_index_files.toList().set { blast_index_list }
process query {
    tag "${prefix}"
    publishDir "${params.output_dir}"

    input:
    set val(sraID), val(sex), val(tissue), file(fasta), file(blast_indexes: "*") from filtered_fastas.combine(blast_index_list)

    script:
    prefix = "${sraID}.${tissue}.${sex}"
    output_blast = "${prefix}.blast.genes.txt"
    """
    echo "threads: \${NSLOTS:-\${threadsMid:-\${NTHREADS:-1}}}"

    blastx -query "${fasta}" -db "${params.blast_index_basename}" -outfmt 6 -max_target_seqs 1 -out "${output_blast}" -num_threads \${NSLOTS:-\${threadsMid:-\${NTHREADS:-1}}}
    """
    // blastx -query "${new}".putative_genes.fa -db uniprot-chicken -outfmt 6 -max_target_seqs 1 -out "${new}".putative_genes.txt -num_threads 10
}


// bedtools intersect -wa -c -a SRR924202.bed -b chickspress_peptides.bed > srr924202_new_screen_intersect.txt
// cut -f 1 SRR924485.putative_genes.txt | sort | uniq > SRR924485.putative_genes.txt.uniq.genes
// grep ">" SRR924485.putative_genes.fa | sed 's/>//g' | sort | uniq > SRR924485.putative_genes.fa.genes
// diff SRR924485.putative_genes.txt.uniq.genes SRR924485.putative_genes.fa.genes | grep ">" | sed 's/>//g' | sed 's/ //g' | sort | uniq > SRR924485.diff.genes
// grep -f SRR924485.diff.genes SRR924485.combined.gtf > SRR924485.combined.lincs.gtf
// cut -f 1 SRR924539.putative_genes.txt | sort | uniq > SRR924539.putative_genes.txt.uniq.genes
// grep ">" SRR924539.putative_genes.fa | sed 's/>//g' | sort | uniq > SRR924539.putative_genes.fa.genes
// diff SRR924539.putative_genes.txt.uniq.genes SRR924539.putative_genes.fa.genes | grep ">" | sed 's/>//g' | sed 's/ //g' | sort | uniq > SRR924539.diff.genes
// grep -f SRR924539.diff.genes SRR924539.combined.gtf > SRR924539.combined.lincs.gtf
// cat SRR924485.combined.lincs.gtf SRR924539.combined.lincs.gtf > Adipose_linc.gtf
// gtf2bed < Adipose_linc.gtf > Adipose_linc.bed
