# TUPaC-lncur:

[other suggestions: TU-CuP LiNC, sounds similar to bioinfo utility cufflinks]
'TUCP-LNCer' is a machine learning approach to identify long noncoding RNAs that encode micropeptides.

## Abstract
Long noncoding RNAs (lncRNAs) are a diverse class of RNAs that have a length of at least 200 nucleotides, and although transcribed from DNA, they are not translated into proteins (or lack an open reading frame of >100 amino acids). Several lncRNAs have been proven to be essential for life, including the ribosomal RNA for protein synthesis, the telomerase RNA in protecting chromosomal ends during mitosis and *Xist* in X-chromosome inactivation in mammalian females. Although categorized as non-protein coding, recent studies have demonstrated that small peptides, known as micropeptides, can be concealed within lncRNA transcripts. LncRNAs containing these micropeptides have been described as Transcripts of Uncertain Coding Potential (TUCP). These micropeptides remained largely undetected in earlier computational analyses due to their short size and non-cononical start sites. However, a micropeptide as short as 11 amino acids, named **torsal-less**, was found to be important for leg development in fruit fly. Additionally, functionally relevant micropeptides have been discovered in human, mouse and chicken. These recent findings raise an important question: are lncRNAs just protein-coding genes that defy traditional concepts of a being protein, such as being at least 100 amino acids in length and evolutionary conservation, that were used in early computational analyses? Here, in this work we aim to identify micropeptides from 'omic' datasets using using a Machine Learning (ML) approach. We have extensively profiled the lncRNA transcriptome in 13 matched, chicken male and female tissues for which RNA-Seq data and protein data is publicaly available using two different pipelines - RMTA () and Evolinc (Ref) and RMTA and a modified EvoTuc (). To emphasize on the discovery of novel micropeptides, we trained the ML algorithm on the publically available 'chickpress' dataset (http://geneatlas.arl.arizona.edu/) that contains mass spectrometry (MS) data and we excluded RNAs that map to known proteins in chicken, and trained ML algorithm on the proteins that overalap with the lncRNAs detected from RNASeq dataset.

## Background
Proteins have been long considered to be the ‘workforce’ for biological systems. Due to the variety in chemical properties of the building blocks of proteins, known as amino acids (AA), proteins are capable of performing a diverse range of functions including enzymatic activity, signalling molecules, structural blocks that hold cells and tissues together, and many more. Despite their central role, the information on DNA that codes for proteins usually represents a small fraction of the genome of higher organisms. For example, the protein-coding genes represent only 2% of the human genome of roughly 3 billion letters; the rest was regarded as ‘junk’ and largely neglected. Recent advances in sequencing technologies have revealed that the ‘junk’ DNA, despite being non-protein coding, was, in fact, transcribed, representing majority of the intracellular pool of RNA. A recent survey found that at least 75% of the human genome was transcribed in at least one cell type (ENCODE Project Consortium 2012). Although some investigators argued that this pervasive transcription could be simply transcriptional noise, classical biochemical studies have firmly established biological roles of a small number of noncoding RNAs, including ribosomal and transfer RNA in protein synthesis, telomerase RNA in protecting chromosomal ends during mitosis and Xist in X-chromosome inactivation in mammalian females (Rinn and Chang 2012). Thus, a broader role of noncoding RNAs in physiology is an intense area of research. 

Noncoding RNAs are broadly classified into two categories based on their size: small (<200 nucleotides) and long noncoding RNAs (lncRNAs; >200 nucleotides). While not being fully translated, some lncRNAs do contain short open reading frames that are translated into micropeptides. In fact, a recent study identified a micropeptide of 46 AA that originated from a lncRNA, was evolutionary conserved between human and mouse and found to regulate muscle physiology via calcium signalling (Anderson et al. 2015). Besides human and mouse, micropeptides are also detected in fruit fly (Pueyo and Couso 2008) and chicken (Cai et al. 2017). A revised analyses of 25 independent studies identified roughly 3,500 TUCPs in human and some of these were indeed translated and detected in proteomic studies (Iyer et al. 2015). 

## Why should we solve it?
Identification of lncRNAs that produce micropeptides will have several advantages. First, the human genome is estimated to house roughly 55,000 lncRNAs, which greatly exceeds the number of protein-coding genes of around 20,000 (Iyer et al. 2015). Despite the large number of lncRNAs discovered rapidly, there is currently insufficient data to comprehend their biological role, devise a classification scheme and  investigate their functions systematically. A subset of lncRNAs that could partially work as micropeptides will immediately allow dissection of their functions using the existing toolset of molecular biology. Furthermore, a list of legitimate micropeptides, as validated by functional assays, might prove to be essential for an improved understanding of the living systems as well as treatment of diseases. The clinical potential of novel micropeptides is immense, ranging from development of novel diagnostic tests, drugs or even vaccines to precision medicine for effective management of patients. 

# Optional things
- Please cite our work -- here is the ICMJE Standard Citation:
- ...and a link to the DOI:
- Awesome Logo

# How to use this software
_*To be filled by computer wizards*_

# Methods
## Software Workflow Diagram
_*Following might have changed*_

![''](https://github.com/NCBI-Hackathons/ncRNA_ML_Features/blob/master/flowchart/jay_flow1.png)


## File structure diagram
#### _Define paths, variable names, etc_

## Installation options:

We provide two options for installing <this software>: Docker or directly from Github.

### Docker

The Docker image contains <this software> as well as a webserver and FTP server in case you want to deploy the FTP server. It does also contain a web server for testing the <this software> main website (but should only be used for debug purposes).

1. `docker pull ncbihackathons/<this software>` command to pull the image from the DockerHub
2. `docker run ncbihackathons/<this software>` Run the docker image from the master shell script
3. Edit the configuration files as below

### Installing <this software> from Github

1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
2. Edit the configuration files as below
3. `sh server/<this software>.sh` to test
4. Add cron job as required (to execute <this software>.sh script)

### Configuration and dependencies

```Examples here```

## Sequence featurization

Non-overlapping subsequences of 100 nt covering the entire length of a lncRNA were featurized using a set of *k*’s to generate *k*-mers (*k* = 2, 3 and 4). For each subsequence, the *k*-mer usage frequencies for each value of *k* was computed, but  the unused *k*-mers, i.e. with frequency = 0, were also included in the vector to ensure equal length outputs. Finally, the *k*-mer usage frequency vector was concatenated into a single vector representing the subsequences.

## Testing

We tested four different tools with <this software>. They can be found in [server/tools/](server/tools/) .

# Additional Functionality

---------




### DockerFile

<this software> comes with a Dockerfile which can be used to build the Docker image.

  1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `cd server`
  3. `docker build --rm -t <this software>/<this software> .`
  4. `docker run -t -i <this software>/<this software>`

### Website

There is also a Docker image for hosting the main website. This should only be used for debug purposes.

  1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `cd Website`
  3. `docker build --rm -t <this software>/website .`
  4. `docker run -t -i <this software>/website`

# References & Resources
(in alphabetical order; see gSheets for sorting etc, then paste in the end)
https://docs.google.com/spreadsheets/d/1aqINk-C1tgVPTUMGZFKQanBXSp_GVIyM7hxrdhtBeMI/edit?usp=sharing


- [Anderson et al. 2015](http://science.sciencemag.org/content/351/6270/271)
- [chickpress dataset on the University of Arizona, USA](http://geneatlas.arl.arizona.edu/)
- [CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features](https://academic.oup.com/nar/article/45/W1/W12/3831091)
- [ENCODE Project Consortium 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3439153/)
- [Evolinc: A Tool for the Identification and Evolutionary Comparison of Long Intergenic Non-coding RNAs](https://www.frontiersin.org/articles/10.3389/fgene.2017.00052/full)
- [Gallus gallus reference proteome](http://www.uniprot.org/proteomes/UP000000539)
- [Iyer et al. 2015](https://doi.org/10.1038/ng.3192)
- [Pueyo and Couso 2008](https://www.sciencedirect.com/science/article/pii/S0012160608011597?via%3Dihub)
- [Rinn and Chang 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3858397/)

