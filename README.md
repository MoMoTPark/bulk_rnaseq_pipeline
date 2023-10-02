### Bulk RNA-seq data processing pipeline with Snakemake
Included `snakefile` in this repository is a general purpose bulk RNA-seq data processing pipeline. I have developed this pipeline with two mapping steps to maximise the sensitivity of finding novel exon junctions (i.e., transcripts) in transcriptomic data. This is a high performance and efficient pipeline designed to be used with AWS EC2 machines that utilises up to 32 cores and 256GB of RAM. The pipeline can be used to process very large datasets with many samples safely without breaking if number of threads and allocated memory is managed properly.