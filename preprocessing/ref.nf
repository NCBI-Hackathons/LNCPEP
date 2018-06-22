process download_ref_fasta {
    storeDir "${params.ref_dir}"

    input:
    val(foo) from Channel.from("foo")

    output:
    file("Gallus_gallus-5.0_genomic.fasta")

    script:
    """
    wget https://de.cyverse.org/dl/d/1BC80779-B57B-4C2F-9557-36D0F934ACB8/Gallus_gallus-5.0_genomic.fasta
    """
    // wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.4_Gallus_gallus-5.0/GCF_000002315.4_Gallus_gallus-5.0_genomic.fna.gz
    // gunzip GCF_000002315.4_Gallus_gallus-5.0_genomic.fna.gz
}

process download_ref_gtf {
    storeDir "${params.ref_dir}"

    input:
    val(foo) from Channel.from("foo")

    output:
    file("Gallus_gallus-5.0_genomic_cleaned.gff")

    """
    wget https://de.cyverse.org/dl/d/9BD182ED-7201-4485-B627-6C6D94192C35/Gallus_gallus-5.0_genomic_cleaned.gff
    """
}
