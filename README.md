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

## Optional pre-requisites

For alignment to a reference genome, you will need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](https://github.com/samtools/samtools) installed. If you only want to annotate/segment/filter the MAS-ISO-seq reads, you may skip this part.

## A fully worked example

```sh
# Download the reads, PacBio index file, the SIRV reference sequence, and the SIRV gene annotations
$ wget gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas10.reads.bam
$ wget gsapubftp-anonymous@ftp.broadinstitute.org:/MasSeqNatBiotech2021/SIRV_MAS_15-10x_mas10.reads.bam.pbi
$ wget https://raw.githubusercontent.com/broadinstitute/MAS-ISO-seq/main/SIRV_Library.fasta
$ wget https://raw.githubusercontent.com/broadinstitute/MAS-ISO-seq/main/SIRV_Library.gff3

# Annotate a chunk of the reads with the appropriate array model (in this case, 'mas10')
$ longbow annotate -m mas10 \
                   -o SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.bam \
                   -c 1/300 \
                   SIRV_MAS_15-10x_mas10.reads.bam
[INFO 2021-10-20 11:11:22 annotate] Invoked via: longbow annotate -m mas10 -o SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.bam -c 1/300 SIRV_MAS_15-10x_mas10.reads.bam
[INFO 2021-10-20 11:11:22 annotate] Running with 11 worker subprocess(es)
[INFO 2021-10-20 11:11:22 annotate] Using The MAS-seq 10 array element model.
[INFO 2021-10-20 11:11:31 annotate] Annotating 2683 reads from chunk 1/300
Progress: 100%|███████████████████████████████████████████████████████████████████████████████████████| 2683/2683 [03:20<00:00, 13.41 read/s]
[INFO 2021-10-20 11:14:51 annotate] Annotated 2683 reads with 107712 total sections.
[INFO 2021-10-20 11:14:51 annotate] Done. Elapsed time: 208.69s. Overall processing rate: 12.86 reads/s.

# Segment the reads according to the annotations
$ longbow segment -m mas10 \
                  -o SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.bam \
                  SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.bam
[INFO 2021-10-20 11:16:29  segment] Invoked via: longbow segment -m mas10 -o SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.bam SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.bam
[INFO 2021-10-20 11:16:29  segment] Running with 11 worker subprocess(es)
[INFO 2021-10-20 11:16:29  segment] Using bounded region splitting mode.
Progress: 0 read [00:00, ? read/s][INFO 2021-10-20 11:16:29  segment] Using The MAS-seq 10 array element model.
[INFO 2021-10-20 11:16:30  segment] Using The MAS-seq 10 array element model.
[INFO 2021-10-20 11:16:46  segment] Segmented 2683 reads with 20835 total segments.
[INFO 2021-10-20 11:16:46  segment] MAS-seq gain factor: 7.77x
[INFO 2021-10-20 11:16:46  segment] Done. Elapsed time: 17.63s.

# Filter out reads that do not conform to the array design
$ longbow filter -m mas10 \
                 -o SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.filtered.bam \
                 SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.bam
[INFO 2021-10-20 11:17:40   filter] Invoked via: longbow filter -m mas10 -o SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.filtered.bam SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.bam
[INFO 2021-10-20 11:17:40   filter] Using The MAS-seq 10 array element model.
[INFO 2021-10-20 11:17:41   filter] Writing reads that conform to the model to: SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.filtered.bam_longbow_filter_passed.bam
[INFO 2021-10-20 11:17:41   filter] Writing reads that do not conform to the model to: SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.filtered.bam_longbow_filter_failed.bam
[INFO 2021-10-20 11:17:41   filter] Filtering according to mas10 model ordered key adapters: Q, C, M, I, O, J, B, D, K, H, R
Progress: 20835 read [00:16, 1259.09 read/s]
[INFO 2021-10-20 11:17:58   filter] Done. Elapsed time: 18.01s.
[INFO 2021-10-20 11:17:58   filter] Total Reads Processed: 20835
[INFO 2021-10-20 11:17:58   filter] # Reads Passing Model Filter: 18350 (88.07%)
[INFO 2021-10-20 11:17:58   filter] # Reads Failing Model Filter: 2485 (11.93%)
[INFO 2021-10-20 11:17:58   filter] Total # correctly ordered key adapters in passing reads: 18350
[INFO 2021-10-20 11:17:58   filter] Total # correctly ordered key adapters in failing reads: 0
[INFO 2021-10-20 11:17:58   filter] Avg # correctly ordered key adapters per passing read: 1.0000 [11]
[INFO 2021-10-20 11:17:58   filter] Avg # correctly ordered key adapters per failing read: 0.0000 [11]

# Align the reads to the SIRV reference sequence
$ samtools fastq SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.filtered.bam_longbow_filter_passed.bam | \
    minimap2 -ayYL -x splice:hq -R "@RG\tID:SIRVmas10\tSM:SIRV" SIRV_Library.fasta - | \
    samtools sort - \
    > SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.filtered.bam_longbow_filter_passed.aligned.bam
$ samtools index SIRV_MAS_15-10x_mas10.reads.annotated.chunk1.segmented.filtered.bam_longbow_filter_passed.aligned.bam
```

![IGV screenshot](https://github.com/broadinstitute/MAS-ISO-seq/blob/main/IGV.png?raw=true)
