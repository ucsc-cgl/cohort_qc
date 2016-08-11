# cohort_qc
Quality control of fastq and bam files at cohort level

## Introduction 

This repo contains two scripts. The python script called Dividing_file_im.py which takes 
as input a directory that contains a batch with samples of Genome sequence data. The python 
script follow two steps. First it takes the zip file of the sample and unzips it to later 
take the txt that contains all the data and divide it according to its content. 
The R script named Graph.R takes as input the unzipped divided by content files and imports them to R. 
Then the R script complete a series of test and outputs a report with all the bad samples 
within the batch.

Prerequisites:
1) Must have at least 10 sample because of the ci.median in the R script.

##Install
R studio libraries:
library(plyr)
library(ggplot2)
library(tcltk)
library(asbio)
library(nortest)
library(tidyr)
library(gridExtra)

To be able to use tidyr you must download XQuartz (for MacOS)

Python 2.7.11

## Run

First run the python script by writing python Directory/of/python/script Directory/of/file
Example:
python cohort_qc/Dividing_file_im.py /Users/Ilian/Desktop/fastQC_treehouse_results
Second run the R script by writing Rscript --vanilla Graph.R Directory/of/split/files NumberOFSlpitSample
Example:
Rscript --vanilla Graph.R /Users/Ilian/Desktop/fastQC_treehouse_results 30

