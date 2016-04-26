source ~/.bashrc

#for reads in synthetic ../../lambda_raw_reads_all/Amplicon_4070-5989
for reads in ../../../reads/lambda_raw_reads_all/Amplicon_4070-5989 synthetic
do

for norm in n y 
do

for offset in 100 # 50 75 200  # 100 # 100 200 300 400 500 600
do

for winSz in 4 8 10 12 16 32 64 128 256 512 1024 2048
do


py mlpy_dtw_align.py \
	../../../referenceGenomes/lambda.fasta \
	lambda_template.model \
	lambda_complement.model \
	$reads \
	$winSz \
	4070 \
	5989 \
	$offset \
	$norm 

done
done
done
done 





	#lambda-2d-reads \
	#test_reads_Amplicon_1 \
	#../../lambda_raw_reads_all/Amplicon_3/ \
# 1 52 1980
# 2 2065 3965
# 3 4070 5989 
