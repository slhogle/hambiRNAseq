#!/usr/bin/env bash
#SBATCH -J preprocess
#SBATCH --time 0-2:00:00
#SBATCH -p serial
#SBATCH --mem 32GB
#SBATCH -n 8

module load fastqc/0.11.8.lua

MYLIB=$SLURM_ARRAY_TASK_ID

ADAPTERS="/homeappl/home/hoglesha/appl_taito/apps/bbmap-v38.61b/resources/adapters.fa"
ARTEFACTS="/homeappl/home/hoglesha/appl_taito/apps/bbmap-v38.61b/resources/sequencing_artifacts.fa.gz"
PHIX="/homeappl/home/hoglesha/appl_taito/apps/bbmap-v38.61b/resources/phix174_ill.ref.fa.gz"

# Remove optical duplicates only. Since these are transcriptomes there will
# be a fairly large amount of duplication in these samples. Especially since
# rRNA depletion wasn't performed

# these settings optimized for NovaSeq. This RNAseq done with S4 flowcell
clumpify.sh dedupe=t optical=t dist=12000 spantiles=f \
in=rawreads/${MYLIB}_R1.fastq.gz \
in2=rawreads/${MYLIB}_R2.fastq.gz \
out=processedreads/${MYLIB}.00-clumped.fq.gz
cd processedreads; rm -f ${MYLIB}.tmp.fq.gz; ln -s ${MYLIB}.00-clumped.fq.gz ${MYLIB}.tmp.fq.gz; cd ..

# Remove low-quality regions
# dont know if this is necessary... Omitting to keep more data
# filterbytile.sh \
# in=processedreads/${MYLIB}.tmp.fq.gz \
# out=processedreads/${MYLIB}.01-tiled.fq.gz
# cd processedreads; rm -f ${MYLIB}.tmp.fq.gz; ln -s ${MYLIB}.01-tiled.fq.gz ${MYLIB}.tmp.fq.gz; cd ..

# trim adapters
bbduk.sh ow=t ref=adapters ordered \
ktrim=r k=23 mink=11 hdist=1 hdist2=1 tbo tpe ftm=5 ordered \
ref=${ADAPTERS} \
in=processedreads/${MYLIB}.tmp.fq.gz \
out=processedreads/${MYLIB}.02-adaptertrim.fq.gz
cd processedreads; rm -f ${MYLIB}.tmp.fq.gz; ln -s ${MYLIB}.02-adaptertrim.fq.gz ${MYLIB}.tmp.fq.gz; cd ..

# filter matches to phiX spike in and common sequencing artefacts
bbduk.sh ow=t k=27 ordered \
ref=${ARTEFACTS},${PHIX} \
in=processedreads/${MYLIB}.tmp.fq.gz \
out=processedreads/${MYLIB}.03-phixfilter.fq.gz
cd processedreads; rm -f ${MYLIB}.tmp.fq.gz; ln -s ${MYLIB}.03-phixfilter.fq.gz ${MYLIB}.tmp.fq.gz; cd ..

# quality trim
bbduk.sh ow=t qtrim=r trimq=6 minlen=70 maxns=2 minavgquality=10 ordered \
in=processedreads/${MYLIB}.tmp.fq.gz \
out=processedreads/${MYLIB}.04-qualitytrim.fastq.gz
cd processedreads; rm -f ${MYLIB}.tmp.fq.gz; ln -s ${MYLIB}.04-qualitytrim.fastq.gz ${MYLIB}.tmp.fq.gz; cd ..

# filter out reads matching rRNA sequences
bbmap.sh ow=t int=f ambiguous=best semiperfectmode=t \
ref=HAMBI-30.16S-sanger_wgs.dedupe.ffn \
ehist=fastqc_postprocess/${MYLIB}.ehist \
qahist=fastqc_postprocess/${MYLIB}.qahist \
mhist=fastqc_postprocess/${MYLIB}.mhist \
idhist=fastqc_postprocess/${MYLIB}.idhist \
scafstats=fastqc_postprocess/${MYLIB}.scafstats \
statsfile=fastqc_postprocess/${MYLIB}.mapstats \
in=processedreads/${MYLIB}.tmp.fq.gz \
outm=processedreads/${MYLIB}.05-rRNA.fastq.gz \
outu=processedreads/${MYLIB}.05-mRNA.fastq.gz

rm processedreads/${MYLIB}.tmp.fq.gz

# get fastq stats of cleaned reads
fastp -i processedreads/${MYLIB}.05-rRNA.fastq.gz \
-j fastqc_postprocess/${MYLIB}_rRNA.json \
-h fastqc_postprocess/${MYLIB}_rRNA.html \
-R ${MYLIB}_rRNA \
--dont_overwrite -w 8 \
--interleaved_in \
-p -P 20 \
--disable_adapter_trimming \
--disable_trim_poly_g \
--disable_quality_filtering \
--disable_length_filtering

fastp -i processedreads/${MYLIB}.05-mRNA.fastq.gz \
-j fastqc_postprocess/${MYLIB}_mRNA.json \
-h fastqc_postprocess/${MYLIB}_mRNA.html \
-R ${MYLIB}_mRNA \
--dont_overwrite -w 8 \
--interleaved_in \
-p -P 20 \
--disable_adapter_trimming \
--disable_trim_poly_g \
--disable_quality_filtering \
--disable_length_filtering

# # get fastq stats of rRNA removed/cleaned reads
# cd fastqc_postprocess
# fastqc -f fastq -o . ../processedreads/${MYLIB}.05-rRNA.fastq.gz
# unzip ${MYLIB}.05-rRNA.zip
# rename "${MYLIB}.05-rRNA_fastqc" "${MYLIB}_rRNA" *
# rm -rf ${MYLIB}_rRNA.zip
# rm -rf ${MYLIB}_rRNA/Icons
# rm -rf ${MYLIB}_rRNA/Images
# rm -rf ${MYLIB}_rRNA/fastqc_report.html
#
# fastqc -f fastq -o . ../processedreads/${MYLIB}.05-mRNA.fastq.gz
# unzip ${MYLIB}.05-mRNA.zip
# rename "${MYLIB}.05-mRNA_fastqc" "${MYLIB}_mRNA" *
# rm -rf ${MYLIB}_mRNA.zip
# rm -rf ${MYLIB}_mRNA/Icons
# rm -rf ${MYLIB}_mRNA/Images
# rm -rf ${MYLIB}_mRNA/fastqc_report.html
# cd ..
