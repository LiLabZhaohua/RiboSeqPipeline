# RiboSeqPipeline
## Scripts for Ribo-seq data analysis by Li Lab
This pipeline was modified from piPipes small RNAseq pipeline: https://github.com/bowhan/piPipes/blob/master/bin/piPipes_smallRNA.sh

To use the pipeline, please install piPipes (https://github.com/bowhan/piPipes) first, then put piPipes_with_rfp file into piPipes main directory (at the same directory with piPipes file).
Before running the pipeline, please follow the instructions of piPipes to install the genomes (for example, by running piPipes install -g GENOME).
Here is an example to call the pipeline, using mouse mm10 genome and 8 threads:
piPipes_with_rfp ribo -i Data.fastq.gz -g mm10 -c 8
