# MAS-ISO-seq

## Increased PacBio transcriptome throughput via programmatic concatenation and statistical deconcatenation

This repository contains links to example data and tools/instructions for initial processing of the data (segmentation, filtering, and alignment).

## Example Data

Much of our validation of MAS-ISO-seq was done using Spike-in RNA Variant Control Mix sequences ([SIRVs set 4](https://www.lexogen.com/sirvs/sirv-sets/) from [Lexogen](https://www.lexogen.com/)).  Two example SIRV datasets (with two files each) are available for download by anonymous FTP:

| Array design | Dataset          | Size  | Download command                                                                                           |
|--------------|------------------|-------|------------------------------------------------------------------------------------------------------------|
| mas10        | SIRV set 4 reads | 17.0G | wget gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas10.reads.bam      |
| mas10        | SIRV set 4 index | 7.0M  | wget gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas10.reads.bam.pbi  |
| mas15        | SIRV set 4 reads | 98.8G | wget gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas15.reads.bam      |
| mas15        | SIRV set 4 index | 31.8M | wget gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas15.reads.bam.pbi  |

\* - The SIRV samples were prepared with two library preparation techniques: a length 10 MAS-ISO-seq array and a length 15 MAS-ISO-seq array.  They were multiplexed into a single pooled sample and sequenced in a single run on a PacBio Sequel IIe.  Our software package, [Longbow](https://github.com/broadinstitute/longbow/releases/tag/v0.2.2), was then used to demultiplex the single SIRV multiplexed sample into two outputs - one for the length 15 array and one for the length 10 array.

## Installing the MAS-seq processing tool, Longbow

From Github source:

```sh
$ git clone git@github.com:broadinstitute/longbow.git
$ cd longbow
$ python -mvenv venv \
      && . venv/bin/activate \
      && pip install -r dev-requirements.txt \
      && pip install -e .
$ longbow version
0.4.3
```

## Processing 

```sh
$ longbow annotate -m mas10 -o SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.bam -c 1/300 SIRV_MAS_15-10x_mas10.reads.bam
$ longbow segment -m mas10 -o SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.bam SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.bam
$ longbow filter -m mas10 -o SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.filtered.bam SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.bam
```

