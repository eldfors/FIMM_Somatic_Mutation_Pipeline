# FIMM_Somatic_Mutation_Pipeline
This repository hosts the FIMM Somatic Mutation Analysis Pipeline. It identifies and annotates somatic mutations, and outputs detailed mutation data. Designed to assist in cancer genomics research by providing comprehensive insights into tumor genetics. The pipeline is implemented in Python.

## Overview

The pipeline processes a pair of tumor and normal genome sequencing reads pre-alinged to a reference genome, and performs a variant calling analysis to identify somatic mutations. It applies a combination of bioinformatics tools including VarScan, SnpEff, and Samtools.

## Dependencies

The FIMM Somatic Mutation Pipeline requires the following tools to be installed and configured:

1. [Python 2.7](https://www.python.org/download/releases/2.7/)
    - Additional Python modules:
      - pandas
      - PyVCF
      - configparser

2. [VarScan](http://dkoboldt.github.io/varscan): VarScan is a platform-independent mutation caller for targeted, exome, and whole-genome resequencing data generated on Illumina, and similar instruments.

3. [SnpEff](http://snpeff.sourceforge.net/): SnpEff is a genetic variant annotation and effect prediction toolbox. It annotates and predicts the effects of variants on genes (such as amino acid changes).

4. [Samtools](http://www.htslib.org/): Samtools is a suite of programs for interacting with high-throughput sequencing data. It allows for efficient manipulation of alignments in the SAM/BAM format, including sorting, indexing, and generating alignments in a per-position format.

4. Reference Genome: The pipeline uses a reference genome for alignment and variant calling. This should be downloaded from [ENSEMBL](http://www.ensembl.org/info/data/ftp/index.html). The specific version (e.g., hg38 for human) will depend on your research needs.

Please follow the respective links to access the official documentation and download pages for each tool.

## Installation and Setup
To install the pipeline, clone this repository to your local machine using git:

First, clone this repository to your local machine using git clone. Replace [url] with the URL of your repository:
git clone https://github.com/eldfors/FIMM_Somatic_Mutation_Pipeline.git

    ...
  
## MIT License

Copyright (c) 2023 Samuli Eldfors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
