Files associated with CMH and file editing post CMH for downstream FLAM.

cmh_metabolomics.R
R code for detecting statistical differences between three treatments of 10 replicated populations in SNP frequency. The input for this code was a snp_table generated from aligned bam files. The snp_table is too large to upload to Github.

meta_20cov_sigSNP.txt
Text file of the significantly differentiated SNPs and their respective positions and frequencies. 

Obtain_Low_PValue_Position_Annotated.py
Python script to extract the lowest p-value SNP within a 50kb region to use in the downstream FLAM analysis.

outputmeta_20cov_sigSNP.txt
Text file of the lowest p-value per 50kb region that contains at least 3 significant SNPs produced by "Obtain_Low_PValue_Position_Annotated.py.
