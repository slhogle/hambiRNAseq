#!/usr/bin/bash

# run all steps in the analysis

# Must untar rawData/16SAmplicon/bbmapRPKM first

# SBW25 variant analysis

# LGPR

# Amplicon analysis
#Rscript ../R/rpkm2tab.R
Rscript ../R/makePhyloseq.R
Rscript ../R/makeFigS2.R
Rscript ../R/shannonDiversity.R
Rscript ../R/ordination.R
#Rscript ../R/corncobEvolution.R 
#Rscript ../R/corncobPredation.R
#Rscript ../R/corncobInteraction.R
Rscript ../R/corncobSBW25.R
Rscript ../R/Rscript makeFig3.R