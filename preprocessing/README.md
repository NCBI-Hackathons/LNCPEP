# Preprocessing Workflow

Data workflow for generating lincRNA's to use for training

# Installation

Install Nextflow in the local directory (requires Java 8)

```
make install
```

Download reference data

```
make ref
```

# Run

```
make run
```

# Reqiurements

- reference genome .fasta

- reference genome annotations .gff

- reference proteome .fasta

- SRA identifier for research sample(s)

# Workflow

Processing steps that will be performed by the Nextflow pipeline

## Setup Steps

- create HiSat2 indexes for reference genome fasta

- create BLAST databases & indexes for proteome fasta

- align SRA reads with `hisat2`

- sort aligned reads with `sambamba`

- assemble transcripts with `stringtie`

- filter low coverage transcripts

- extract fasta sequences for transcripts with `gffread`

- filter fasta transcripts for sequences >200bp

- query transcript fasta sequences against proteome database with `blastx` to identify protein coding transcripts

- 
