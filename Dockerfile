FROM jupyter/datascience-notebook:r-4.0.3
MAINTAINER Edwin Juarez <ejuarez@ucsd.edu>

ENV LANG=C LC_ALL=C
USER root

# I preffer to have all the code here rather than in an Rscript file:
# Remember that you don't need the ; after the if statement.
RUN Rscript -e  'install.packages("foreign", repos = "https://cloud.r-project.org/", quiet = TRUE);\
install.packages("Hmisc", repos = "https://cloud.r-project.org/", quiet = TRUE);\
install.packages("getopt", repos = "https://cloud.r-project.org/", quiet = TRUE);\
install.packages("optparse", repos = "https://cloud.r-project.org/", quiet = TRUE);\
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org/", quiet = TRUE);\
BiocManager::install("tximport", quiet = TRUE);\
BiocManager::install("DESeq2", quiet = TRUE);\
BiocManager::install("GenomicFeatures", quiet = TRUE);\
BiocManager::install("rhdf5");\
sessionInfo()'

# build using this:
# docker build -t genepattern/tximport_deseq2:1.0 .
