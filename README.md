# RiboSeqPipeline
## Scripts for Ribo-seq data analysis by Li Lab
This pipeline was modified from piPipes small RNAseq pipeline: https://github.com/bowhan/piPipes/blob/master/bin/piPipes_smallRNA.sh

To use the pipeline, please install piPipes (https://github.com/bowhan/piPipes) first, then put piPipes_with_rfp file into piPipes main directory (at the same directory with piPipes file).

Before running the pipeline, please follow the instructions of piPipes to install the genomes (for example, by running piPipes install -g GENOME).

Here is an example to call the pipeline, using mouse mm10 genome and 8 threads:

piPipes_with_rfp ribo -i Data.fastq.gz -g mm10 -c 8

This pipeline has beed used in our recent publication:

Sun, Y.H., Zhu, J., Xie, L.H. et al. Ribosomes guide pachytene piRNA formation on long intergenic piRNA precursors. Nat Cell Biol 22, 200–212 (2020). https://doi.org/10.1038/s41556-019-0457-4
