personal script depository for RNA-seq analysis \
Zihan Wang, 2025 
## Raw data Processing
**Input**: config.yaml, Snakefile, ref genome(fa or fa.gz), annotation (gff, gtf, gff3), may need hand-input strandness \
--------------------------------- \
	trimming: trimgalore \
	post trimming QC: FastQC \
	alignment：STAR \
----------------------------------\
**Output**:(1)featureCount.txt (2)RPKM bigwig (3)MultiQC 

## DEG ANAlysis
*Standard DEG*: in progress \
*Linear Interaction Model*: in progress 
