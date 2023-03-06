#!/bin/bash

# If you just want to follow along without really building a whole genome
# graph, from the vg test directory, you can set:
CONSTRUCTION=./wiki-data/mock-hs37d5
LARGE_DISK_TMP="$(dirname "$(mktemp -u)")"

DATA=${CONSTRUCTION}/data
REFERENCE=${DATA}/hs37d5.fa
GRAPHS=${CONSTRUCTION}/graphs
MAPPING=${GRAPHS}/node_mapping

mkdir -p ${GRAPHS}

LOGS=${CONSTRUCTION}/logs
mkdir -p ${LOGS}
LOGFILE=${LOGS}/vg_all.log
rm -f $LOGFILE

START=$(date +%s)
echo "Start time: $START" >> $LOGFILE

# Build the VG graphs
(seq 1 22; echo X; echo Y) | parallel -j 24 \
    "vg construct -r $REFERENCE -v ${DATA}/chr{}.vcf.gz -R {} -C -a > ${GRAPHS}/chr{}.vg"

STOP=$(date +%s)
echo >> $LOGFILE
echo "Checkpoint: $STOP" >> $LOGFILE
echo "Construction: $(($STOP-$START)) seconds" >> $LOGFILE

# Harmonize the node ids
vg ids -j -m $MAPPING $(for i in $(seq 1 22; echo X; echo Y); do echo ${GRAPHS}/chr${i}.vg; done)

STOP=$(date +%s)
echo >> $LOGFILE
echo "Finish time: $STOP" >> $LOGFILE
echo "Total time: $(($STOP-$START)) seconds" >> $LOGFILE

while read l; do bsub -J $l -e ${l}.log -o ${l}.log -n 4 -q high "vg giraffe --gbz-name auto_ref_vcf.giraffe.gbz --minimizer-name auto_ref_vcf.min --dist-name auto_ref_vcf.dist --progress --fastq-in reads/${l}.rd1.fq.gz --fastq-in reads/${l}.rd2.fq.gz --sample ${l} --max-multimaps 5 -t 4 > ${l}_M5.gam"; done < reads/samples.list

vg giraffe --gbz-name chr10_1st20M_allsv.gbz --progress --fastq-in /public/home/stgui/work/00_TEO_Genomes/07_main/04_SV_Geno_EVO/00_akvw_sv/06_vg/99_benchmark/01_sim_20X_read_chr10/Zm-MO17.rd1.fq.gz --fastq-in /public/home/stgui/work/00_TEO_Genomes/07_main/04_SV_Geno_EVO/00_akvw_sv/06_vg/99_benchmark/01_sim_20X_read_chr10/Zm-MO17.rd2.fq.gz --sample Zm-MO17 --max-multimaps 5 -t 4 > Zm-MO17/Zm-MO17_M5_inter.gam