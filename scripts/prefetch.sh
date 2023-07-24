#!/bin/bash
#$0 for extra params: like token, --ngc meta_data/prj.ngc
export PATH=/Share2/home/lulab1/APP/sratoolkit.2.11.2-centos_linux64/bin/:$PATH

## check if token file exist (needed for EGA and dbGaP)
#need newest version (2.30.0)
if [ -e ./meta_data/*.ngc ];then
	ngc_file=`ls ./meta_data/*.ngc | head -n1`
	ngc="--ngc ${ngc_file}"
else
	ngc=" "
fi

## 1.prefetch
while [ `ls ./sra/*.sra | wc -l` != `cat ./meta_data/sample_ids.txt |wc -l` -a `ls ./sra/*/*.sra | wc -l` != `cat ./meta_data/sample_ids.txt |wc -l` ]; 
do
prefetch -v --max-size 500G ${ngc} --option-file ./meta_data/sample_ids.txt -O ./sra;
done;

## 2.dump 
#note: 
#1.that the outfiles may not always be fastq files, and the file number may change
#2.sra file may be in a sample dir inside the ./sra dir, like ./sra/SRR00009/SRR00009.sra
sample1=`head -n1 ./meta_data/sample_ids.txt`
if [ -e ./sra/*${sample1}*.sra ];then
	echo "./sra/*.sra exist, will dump to ./fastq"

	for i in  `cat ./meta_data/sample_ids.txt`
	do echo $i ;done | parallel -k -I {} -j 3 " \
        echo \"{} start at `date`!\"; \
	fastq-dump -v --split-3 --gzip ./sra/{}.sra -O ./fastq/; \
	echo \"{} done at `date`!\" "
	
	echo "list: `cat ./meta_data/sample_ids.txt|wc -l`";
	echo "sra: `ls ./sra|wc -l`";
	echo "fastq: `ls ./fastq|wc -l`"
elif [ -e ./sra/$sample1/*.sra ];then
        echo "./sra/${sample1}/*.sra exist, will dump to ./fastq"

        for i in  `cat ./meta_data/sample_ids.txt`
        do echo $i ;done | parallel -k -I {} -j 3 " \
        echo \"{} start at `date`!\"; \
        fastq-dump -v --split-3 --gzip ./sra/{}/{}.sra -O ./fastq/; \
        echo \"{} done at `date`!\" "

        echo "list: `cat ./meta_data/sample_ids.txt|wc -l`";
        echo "sra: `ls ./sra|wc -l`";
        echo "fastq: `ls ./fastq|wc -l`"
else
	echo "./sra/*.sra not exist, will not dump"
fi
