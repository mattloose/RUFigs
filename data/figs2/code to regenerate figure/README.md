# "Supplementary Figure 2" of the paper.
# Scripts and data files to generate the plots.
# Run the following .. 


#!/bin/bash

# Edit the following paths as need be...
pyPath="/cygdrive/c/Python27"
rPath="/cygdrive/c/Program\ Files/R/R-3.2.2/bin/x64"


# Generate model files from one of the read files ....
${pyPath}/python getmodels.py --read Amplicon_4070-5989/llssbzms2p35x_20151004_readuntiludududududu_RU21_lambdaPCR_2922_1_ch72_file25_strand.fast5


# Run the Analysis ....
sh run-lambda-2d.sh ${pyPath} | tee out


# Generate the Figures ...
${rPath}/R.exe --no-save --vanilla < run.r
