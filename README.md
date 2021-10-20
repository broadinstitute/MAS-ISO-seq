# MAS-ISO-seq

## Increased PacBio transcriptome throughput via programmatic concatenation and statistical deconcatenation

This repository contains links to example data and tools/instructions for initial processing of the data (segmentation, filtering, and alignment).

## Example Data

Two example datasets (with two files each) are available for download by anonymous FTP:

There were two datasets from this study: 

| Array design | Dataset          | size  | URL                                                                                                   |
|--------------|------------------|-------|-------------------------------------------------------------------------------------------------------|
| mas10        | SIRV set 4 reads | 17.0G | gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas10.reads.bam      |
| mas10        | SIRV set 4 index | 7.0M  | gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas10.reads.bam.pbi  |
| mas15        | SIRV set 4 reads | 98.8G | gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas15.reads.bam      |
| mas15        | SIRV set 4 index | 31.8M | gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas15.reads.bam.pbi  |

\* - The SIRV samples were prepared with two library preparation techniques: a length 10 MAS-ISO-seq array and a length 15 MAS-ISO-seq array.  They were multiplexed into a single pooled sample and sequenced in a single run on a PacBio Sequel IIe.  Our software package, [Longbow](https://github.com/broadinstitute/longbow/releases/tag/v0.2.2), was then used to demultiplex the single SIRV multiplexed sample into two outputs - one for the length 15 array and one for the length 10 array.

## Installing the MAS-seq processing tool, Longbow

From Github source:

```sh
$ git clone git@github.com:broadinstitute/longbow.git \
      && cd longbow
Cloning into 'longbow'...
remote: Enumerating objects: 999, done.
remote: Counting objects: 100% (585/585), done.
remote: Compressing objects: 100% (330/330), done.
remote: Total 999 (delta 365), reused 418 (delta 255), pack-reused 414
Receiving objects: 100% (999/999), 5.52 MiB | 10.52 MiB/s, done.
Resolving deltas: 100% (508/508), done.

$ python -mvenv venv \
      && . venv/bin/activate \
      && pip install -r dev-requirements.txt \
      && pip install -e .
Collecting isort==4.3.21 (from -r dev-requirements.txt (line 1))
  Using cached https://files.pythonhosted.org/packages/e5/b0/c121fd1fa3419ea9bfd55c7f9c4fedfec5143208d8c7ad3ce3db6c623c21/isort-4.3.21-py2.py3-none-any.whl
Collecting tox (from -r dev-requirements.txt (line 2))
  Using cached https://files.pythonhosted.org/packages/78/b0/ce98616ec9c3f270495a2493cde4d81b1f499057222ae77a8103aea59777/tox-3.24.4-py2.py3-none-any.whl
Collecting wheel (from -r dev-requirements.txt (line 3))
  Using cached https://files.pythonhosted.org/packages/04/80/cad93b40262f5d09f6de82adbee452fd43cdff60830b56a74c5930f7e277/wheel-0.37.0-py2.py3-none-any.whl
(...snip...)
      
$ longbow version
0.4.3
```
