# TUCPa LINCer: 
'TUCPa LINCer' is a machine learning approach to identify  long noncoding RNAs that has a potential to encode micropeptides. 

Long noncoding RNAs (lncRNAs) are RNAs of at least 200 nucleotides in length produced from the parts of the DNA that traditionally do not make proteins. While not being fully translated, some lncRNAs do contain short open reading frames and are translated into short peptides known as micropeptides. Due to their short size, micropeptides remained largely undetected in earlier computational analyses. However, a micropeptide as short as 11 amino acids, named torsal-less, was found was found to impact leg development in fruit fly. Additionally, functionally relevant micropeptides have been discovered in human, mouse and chicken. These recent findings raise an important question: are lncRNAs just protein-coding genes that defy traditional concepts of a being protein, such as being at least 100 amino acids in length and evolutionary conservation, that were used in early computational analyses? 

Here, we aim to identify micropeptides from 'omic' datasets using using machine learning (ML). We trained the ML algorithms on the publically available 'chickpress' dataset (http://geneatlas.arl.arizona.edu/) that contains RNA-seq and mass spectrometry data from 13 matched, chicken male and female tissues. To emphasize on discovery of novel micropeptides, we excluded RNA that map to known proteins in chicken, and trained ML on the remaining proteins that were predicted to be less than **XX** amino acids in size to identity features of the cognate RNA sequence that might result into protein synthesis. 

# Optional things
- Please cite our work -- here is the ICMJE Standard Citation:
- ...and a link to the DOI:
- Awesome Logo



Overview Diagram

# How to use this software

# Software Workflow Diagram

# File structure diagram
#### _Define paths, variable names, etc_

# Installation options:

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

# Testing

We tested four different tools with <this software>. They can be found in [server/tools/](server/tools/) .

# Additional Functionality

---------

## What's the problem?
Proteins have been long considered as the ‘workforce’ for biological systems. Due to the variety in chemical properties of the building blocks of proteins, known as amino acids, proteins are capable of performing a diverse range of functions such as enzymatic activity, signalling molecules, structural blocks that hold cells and tissues together, and many more. Despite their central role, the information on DNA that codes for proteins usually represent a small fraction of the genome of higher organisms. For example, the protein-coding genes represent only 2% of the human genome ; the rest was regarded as ‘junk’ and largely neglected . 

Recent advances in sequencing technologies have revealed that the ‘junk’ DNA, despite being non-protein coding, was transcribed at high levels, representing majority of the intracellular RNA pool. A recent survey found that at least 75% of the human genome was transcribed in at least one cell type (ENCODE Project Consortium 2012). Although some investigators argued that this pervasive transcription could be simply transcriptional noise, classical biochemical studies have firmly established biological roles of a small number of noncoding RNAs, such as ribosomal and transfer RNA in protein synthesis, telomerase RNA in protecting chromosomal ends during mitosis and Xist in X-chromosome inactivation in mammalian females (Rinn and Chang 2012). Thus, a broader role of noncoding RNAs in physiology is an intense area of research. 

Noncoding RNAs are broadly classified into two categories based on their size: small (<200 nucleotides) and long noncoding RNAs (lncRNAs; >200 nucleotides). While not being fully translated, some lncRNAs do contain short open reading frames that are translated into peptides. In fact, a recent study identified a micropeptide of 46 amino acids (aa) that originated from a lncRNA, was evolutionary conserved between human and mouse and found to regulate muscle physiology via calcium signalling (Anderson et al. 2015). Besides human and mouse, micropeptides are also detected in fruit fly (Pueyo and Couso 2008) and chicken (Cai et al. 2017). 

## Why should we solve it?
Functional role of micropeptides raises an important question: are lncRNAs just protein-coding genes that defy traditional concepts of a being protein, such as being at least 100 aa in length and evolutionary conservation, that were used in early computational analyses? A revised analyses of 25 independent studies identified ~3,500 transcripts of unknown coding potential (TUCPs) in human and some of these were indeed translated and detected in proteomic studies (Iyer et al. 2015). Thus, identifying lncRNAs that are TUCPs 



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

- [Evolinc: A Tool for the Identification and Evolutionary Comparison of Long Intergenic Non-coding RNAs](https://www.frontiersin.org/articles/10.3389/fgene.2017.00052/full)

- [CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features](https://academic.oup.com/nar/article/45/W1/W12/3831091)
