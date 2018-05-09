# A thermoresponsive genetic network mediating control of bud break in hybrid aspen

## Authors
**Rajesh Kumar Singh1, Jay Prakash Maurya1, Abdul Azeez1, 2, Pal Miskolczi1, Szymon Tylewicz1, 3, Nicolas Delhomme1, Katja Stojkovic1, Viktor Busov2, Rishikesh P. Bhalerao1**

## Affiliations

1 Umeå Plant Science Centre, Department of Forest Genetics and Plant Physiology, Swedish University of Agricultural Sciences, SE-901 87 Umeå, Sweden
2 School of Forest Resources and Environmental Science, Michigan Technological University, Houghton, MI 49931
3 Department of Plant and Microbial Biology, University of Zürich, Zollikerstrasse 107, 8008 Zürich, Switzerland

## Abstract

In boreal and temperate ecosystems, temperature signal regulates the reactivation of growth (bud break) in perennials in the spring. Molecular basis of temperature mediated control of bud break is poorly understood. We have identified a genetic network mediating the control of bud break in model tree hybrid aspen. The key components of this network are transcription factors SVL (SHORT VEGETATIVE PHASE-LIKE), closely related to Arabidopsis floral repressor SHORT VEGETATIVE PHASE, and its downstream target TCP18, a tree ortholog of a branching regulator in Arabidopsis. SVL and TCP18 are downregulated by low temperature and genetic evidence, herein presented, demonstrates their role as negative regulators of bud break. SVL mediates bud break by antagonistically acting on gibberellic acid (GA) and Abscisic acid (ABA) pathways, which we demonstrate are positive and negative regulators of bud break, respectively. Thus, our results reveal the mechanistic basis for temperature-cued seasonal control of a key phenological event in perennial plants.

## Content of this repository

This repository contains the custom script (in the src folder) mentioned in the supplementary material of  the manuscript, used for performing the histone methylation analysis of the SVP gene locus.

## Data

The subset of the sequencing data corresponding to the SVL gene locus (+ 1kb upstream/downstream) can be found in the data directory.

## Methods

For the ChIP-seq experiment the apical buds of three biological replicates were collected from WT plants after 10 weeks in SD (10WSD) and after an additional 4 week cold (4WC) treatment. ChIP assays were carried described above. Anti-Trimethyl-Histone H3 (LysK27) (#07-449, Millipore) and Anti-Histone H3 (ab1791, Abcam) antibodies were used for chromatin immunoprecipitation. Ovation Ultralow IL Multiplex System I (Part No. 0304, NuGEN) was used to generate the sequencing libraries according to the product instructions. Pair end sequencing was done by BGI-Tech.

Sequencing reads were processed following our guidelines<sup>1</sup>. Briefly, reads quality was first assessed using FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), v0.11.4. Reads mapping to ribosomal RNA (rRNA) were quantified and filtered using SortMeRNA<sup>2</sup> (v2.1;settings `--log --paired_in --fastx--sam --num_alignments 1`) using the rRNA sequences provided with SortMeRNA. Reads were then filtered to remove adapters and trimmed for quality using Trimmomatic<sup>3</sup> (v0.36 ; settings `TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50`). After every filtering step, FastQC was run again to ensure that no technical artefacts were introduced.  Reads were then mapped to the hybrid aspen genome (Populus tremula x tremuloides, clone T89) using STAR<sup>4</sup> with settings `--outQSconversionAdd -31 --outReadsUnmapped Fastx`. Reads were later on remapped using BWA-MEM <sup>5</sup> with default settings to comparable results. Peaks were called genome wide using MACS2<sup>6</sup> with the following non-default parameters: `-f BAM -g 2.7e8 -s 45 --verbose 3 --nomodel --shiftsize 100 --to-large --keep-dup all`, on sequencing libraries down-sampled to 10 million PE reads. This down-sampled library depth (10M) had been estimated by an ad-hoc saturation/rarefaction analysis based on the number of peaks identified by MASC2 in varied subsets of the original dataset. These downstream analyses (peak-calling, saturation, etc.) were solely used to estimate the fraction of the genome mapped under the different growing conditions. The obtained ratios were used as part of the data normalization for the analysis of the SVL locus histone methylation status.
Reads mapped to the sequence of SVL gene including 1 kb upstream and downstream region were extracted from the alignment. Coverage in the above region was calculated, log2 transformed and corrected for the abundance differences between samples (i.e. the fraction of the genome mapped under the different growing conditions in the 10M PE read subset; the latter selection addressing any library size factor scaling otherwise required). Finally, the H3K27me3 abundance was normalized by H3 abundance. These were used to compare differences in histone methylation between the two time points.

## Reference
<sup>1</sup> Delhomme _et al._ http://www.epigenesys.eu/en/protocols/bio-informatics/1283-guidelines-for-rna-seq-data-analysis
<sup>2</sup> Kopylova _et al._, Bioinformatics, 2012
<sup>3</sup> Bolger _et al._, Bioinformatics, 2014
<sup>4</sup> Dobin _et al._, Bioinformatics ,2013
<sup>5</sup> Heng Li, arXiv, 2013
<sup>6</sup> Zhang _et al._, 2008
