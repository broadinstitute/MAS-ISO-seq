# MAS-ISO-seq

## Increased PacBio transcriptome throughput via programmatic concatenation and statistical deconcatenation

This repository contains links to example data and tools/instructions for initial processing of the data (segmentation, filtering, and alignment).

## Example Data

All data from this study are available online (or are in the process of being uploaded).  

There were two datasets from this study: 

| Array design | Dataset          | size  | URL                                                                                                   |
|--------------|------------------|-------|-------------------------------------------------------------------------------------------------------|
| mas10        | SIRV set 4 reads | 17.0G | [SIRV_MAS_15-10x_mas10.reads.bam](gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas10.reads.bam)      |
| mas10        | SIRV set 4 index | 7.0M  | [SIRV_MAS_15-10x_mas10.reads.bam.pbi](gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas10.reads.bam.pbi) |
| mas15        | SIRV set 4 reads | 98.8G | [SIRV_MAS_15-10x_mas15.reads.bam](gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas15.reads.bam)     |
| mas15        | SIRV set 4 index | 31.8M | [SIRV_MAS_15-10x_mas15.reads.bam.pbi](gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas15.reads.bam.pbi) |

\* - The SIRV samples were prepared with two library preparation techniques: a length 10 MAS-ISO-seq array and a length 15 MAS-ISO-seq array.  They were multiplexed into a single pooled sample and sequenced in a single run on a PacBio Sequel IIe.  Our software package, [Longbow](https://github.com/broadinstitute/longbow/releases/tag/v0.2.2), was then used to demultiplex the single SIRV multiplexed sample into two outputs - one for the length 15 array and one for the length 10 array.  These demultiplexed files are what is currently available in the Terra workspace.
