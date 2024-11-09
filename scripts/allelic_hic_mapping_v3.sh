#!/bin/bash

### Hi-C allelic mapping for Ming Hu's lab
### required: 
### pairtools v1.0.2
### python packages: os, gzip, numpy
### R packages: R.utils, data.table, dplyr

usage() { 
cat <<EOF
Usage: bash allelic_hic_mapping_v3.sh [-h] [-i input.bam] [-o out_dir] [-c chromsize] [-p1 phase1] [-p2 phase2]
Description: allelic mapping of Hi-C reads
Options:    
    -h           Print help and exit
    -i           input bam files 
    -c           chromsize files. Need to be in "chr_allele" format: "chr1_CC021      194003716" 
    -o           output path
    -p1			 allele1?
    -p2			 allele2?
EOF
    exit 1
}

while getopts ":i:c:o:p1:p2:h" flag
do
    case "${flag}" in
        i) f=${OPTARG};;
        o) out_path=${OPTARG};;
		c) chromsize=${OPTARG};;
		p1) phase1=${OPTARG};;
		p2) phase2=${OPTARG};;
		h) usage;;
    esac
done

if [ -z "${f}" ] || [ -z "${out_path}" ] || [ -z "${chromsize}" ]; then
    usage
fi

if [ -z "${phase1}" ] || [ -z "${phase2}" ]; then
    phase1=$(grep 'chr1_' ${chromsize} | awk -F'[_\t]' '{print $2}' | head -1)
    phase2=$(grep 'chr1_' ${chromsize} | awk -F'[_\t]' '{print $2}' | tail -1)
fi

# bwa mem -SP -t 16 ${ref} ${s}_1.fastq.gz ${s}_2.fastq.gz | samtools view -bhS - > ${map_dir}/${s}.bam

echo "Parsing Hi-C pairs..."
s=`basename $f .bam` 
pairtools parse --min-mapq 0 --add-columns XA,NM,AS,XS --nproc-in 16 --nproc-out 16 --drop-sam --walks-policy all -c ${chromsize} -o ${out_path}/${s}.unphase.pairs.gz ${f}
python allelic_mphase_v3.py --in ${out_path}/${s}.unphase.pairs.gz --out ${out_path}/${s}.phase --phase1 ${phase1} --phase2 ${phase2}

echo "Splitting phased pairs..."
for sfile in ${out_path}/${s}.phase*pairs
do
	fname=`basename $sfile .pairs`
	bgzip ${sfile}
	pairtools sort --nproc-in 16 --nproc-out 16 -o ${out_path}/${fname}.sorted.pairs.gz ${sfile}.gz
	pairtools dedup --mark-dups --nproc-in 16 --nproc-out 16 --extra-col-pair phase1 phase2 --output-stats ${out_path}/${fname}.dedup.stats -o ${out_path}/${fname}.dedup.pairs.gz ${out_path}/${fname}.sorted.pairs.gz
	zcat ${out_path}/${fname}.dedup.pairs.gz | grep -v '#' | cut -f17-18 | sort | uniq -c > ${out_path}/${fname}.dedup.pairs.info
	rm ${out_path}/${fname}.sorted.pairs.gz
	Rscript phc.summarize_pairs_lec.R ${out_path}/${fname}.dedup.stats
done
