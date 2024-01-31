
# Summit

![Maintainer](https://img.shields.io/badge/maintainer-gartician-blue)

Summit is a re-interpretation of [GoPeaks](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02707-w) written in Python to call narrow and broad peaks in epigenetic CUT&Tag/CUT&RUN/ATAC-Seq datasets. The major changes in this refactor is a signficant reduction of computational demands (please see below). If you use GoPeaks or Summit in your publication, please cite the original paper: 

Yashar, W.M., Kong, G., VanCampen, J. et al. GoPeaks: histone modification peak calling for CUT&Tag. Genome Biol 23, 144 (2022). https://doi.org/10.1186/s13059-022-02707-w

# Quick-Start

```bash
$ python summit.py -h

usage: summit.py [-h] -b BAM -o OUTPUT [-c CONTROL] [-m MINREADS] [-d MINDIST] [-w MINWIDTH] [-t STEP] [-l SLIDE] [-p PVAL] [-s CHROMSIZE] [--broad] [--chromosomes CHROMOSOMES] [--canonical [CANONICAL]] [-v] [--version]

Summit is a peak caller designed for CUT&TAG/CUT&RUN sequencing data. 
Summit by default works best with narrow peaks such as H3K4me3 and
transcription factors. Summit can be used with the "--broad" flag to
call broad peaks like H3K27Ac/H3K4me1. We encourage users to explore
the parameters of Summit to analyze their data.
        

options:
  -h, --help            show this help message and exit

Required Arguments:
  -b BAM, --bam BAM     Input BAM file with histone mark or protein-of-interest
  -o OUTPUT, --output OUTPUT
                        Output prefix e.g. data/peaks/sample-1

Summit Parameters:
  -c CONTROL, --control CONTROL
                        Control BAM file
  -m MINREADS, --minreads MINREADS
                        Analyze bins greater than mindreads. Default: 15
  -d MINDIST, --mindist MINDIST
                        Merge bins within <mindist> base pairs. Default: 1000
  -w MINWIDTH, --minwidth MINWIDTH
                        Remove peaks less than <mindwidth> base pairs wide. Default: 150
  -t STEP, --step STEP  Bin size for genome bins. Default: 100
  -l SLIDE, --slide SLIDE
                        Slide size for genome bins. Default: 50
  -p PVAL, --pval PVAL  Filter away bins above <pval> significance level. Default: 0.05
  -s CHROMSIZE, --chromsize CHROMSIZE
                        Tab-delimited file containing chromosome name and sizes.

Optional Arguments:
  --broad               Call broad peaks with step=5000, slide=1000, mindist=3000.
  --chromosomes CHROMOSOMES
                        Use a regex string to subset chromosomes.
  --canonical [CANONICAL]
                        Subset canonically-named chromosomes. Default: 'chr[0-9XY]+$'
  -v, --verbose         Run the program in verbose mode
  --version             Display Summit version
```

# Recommended Parameters

| Sequencing Modality | Protein of Interest                             | Recommended Parameters |
|---------------------|-------------------------------------------------|------------------------|
| CUT&Tag or CUT&RUN  | Transcription factor                            | Default parameters     |
| CUT&Tag or CUT&RUN  | Narrow histone mark (H3K4me3)                   | Default parameters     |
| CUT&Tag or CUT&RUN  | Broad histone mark (H3K4me1, H3K27Ac, H3K27me3) | `--broad`              |
| ATAC-Seq            | Open chromatin								    | Default parameters     |

# Major Changes

Summit is meant to be an ultrafast and light adaptation to the original GoPeaks algorithm, while keeping true to the previous findings. The following benchmarks come from the same samples used in the original paper found in GoPeaks. Broadly, these data contained previously published data (K562 samples from Kaya-Okur et. al. 2019) and in-house cell-line data (Kasumi samples from [GSE190793](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190793)). The following benchmarks were obtained from Oregon Health and Science Univeristy's high-performance clusters (HPCs) Exacloud.

## Lower Memory Requirements

![memory-requirements](doc/max-uss.png)

Max USS tracks the maximal memory used in a process. This metric is commonly employed in HPCs to kill any process that takes up much RAM. While CUT&Tag data frequently use 3-8 million paired-end reads (<200MB), GoPeaks requires a prohibitively high amount of memory (20-30GBs) to process the data. By comparison, Summit typically uses less than 2GB of memory with very little deviation between runs (individual points). This makes Summit between 90-95% more efficient than GoPeaks! 

## Faster Run Times

![runtimes](doc/time-diff.png)

Not only is Summit more memory-efficient than GoPeaks, Summit commonly takes 1/2 to 1/3 the time to process the same dataset regardless of input sequencing depth.

## Lower CPU Times

![cpu-times](doc/cpu-time.png)

In addition to reduced overall processing time, Summit requires less CPU time. This is the time that the CPU is actually executing a command, as opposed to waiting for data input/output.

## Memory and Time Requirements by Different Read Depths

Memory Usage                                 | Time Usage
:----------------------------------------------------:|:-----------------------------------------------------:
![K562-memory](doc/all-mem-benchmark.png)  | ![K562-time](doc/all-time-benchmark.png)

If we take an example dataset and downsampled the depth from 90-20%, we can see the effects of memory requirement and program duration over many sequencing depths for a specific (but representative) sample. In the plots above, we can observe the following:

1. Maximal memory consumption is invariant to sequencing depth for Summit. 
2. Summit's program duration increases linearly with sequencing depth. GoPeaks also increases linearly, but with much higher deviation in its trajectory.
3. While the time differences are not significantly different (e.g. 6 minutes vs. 3 minutes), the increased efficiency of Summit over GoPeaks allows users to call peaks on samples with very high sequencing depth, or process hundreds of samples simultaneously without requesting prohibitively high computational resources.

# Notes 

Treated sample read depths ranged between 3M to 30M as show below:

```c
K562_1_H3K27Ac.sorted.markd.bam    4_278_700
K562_1_H3K4me1.sorted.markd.bam    6_867_838
K562_1_H3K4me2.sorted.markd.bam    6_231_210
K562_1_H3K4me3.sorted.markd.bam    7_478_036
K562_1_IgG.sorted.markd.bam        3_073_808
K562_2_H3K27Ac.sorted.markd.bam    5_770_016
K562_2_H3K4me1.sorted.markd.bam    8_550_358
K562_2_H3K4me2.sorted.markd.bam    6_256_098
K562_2_H3K4me3.sorted.markd.bam    7_490_544
K562_2_IgG.sorted.markd.bam        1_891_652
Kasumi_1_H3K27Ac.sorted.markd.bam  11_793_316
Kasumi_1_H3K4me3.sorted.markd.bam  30_640_304
Kasumi_1_IgG.sorted.markd.bam      3_213_274
Kasumi_2_H3K27Ac.sorted.markd.bam  10_169_058
Kasumi_3_H3K27Ac.sorted.markd.bam  5_157_920
Kasumi_3_H3K4me3.sorted.markd.bam  13_668_580
```