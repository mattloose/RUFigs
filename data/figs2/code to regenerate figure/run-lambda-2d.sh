source ~/.bashrc

for reads in synthetic Amplicon_4070-5989
do

for norm in n y 
do

for offset in 50 # 100 # 100 200 300 400 500 600
do

for winSz in 4 8 10 12 16 32 64 128 256 512 1024 2048
do


${1}/python mlpy_dtw_align.py \
	lambda.fasta \
	template.model \
	complement.model \
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
