#!/bin/bash

set -e
set -o pipefail
set -u



if [[ ${#@} != 1 ]]; then
    echo "ERROR: Provide one required argument: configuration file" >&2
    exit 1
fi

CONF=$1

if [[ ! -s ${CONF} ]]; then
    echo "ERROR: Configuration file not found" >&2
    exit 1
fi


###### Variables ######
# Do not change:
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
WD="$(pwd)"

source $CONF
export DOMAIN
export OUTPUT


###### Functions ######
function prefix_anno {
    local prefix=$1
    local fields=$2
    echo ${fields} | awk -v RS=, -v ORS=, '{print "'${prefix}'_"$1":="$1}' | sed 's/,$//'  # "DST_TAG:=SRC_TAG"
}

function select_transcript {
    local chrom=$1
    local filename=$2
    local transcripts
    transcripts=$(cat ${OUTPUT}/transcripts_${DOMAIN}_chr${chrom}.txt | tr \\n , | sed 's/,$//')
    bcftools +split-vep -r ${chrom} -a ANN -p VEP_ -c - -d ${filename} | \
	bcftools annotate -i 'VEP_Feature="'${transcripts}'"' -x INFO/ANN -Ob -o ${filename/.bcf.gz/_transcript.bcf.gz}
    bcftools index ${filename/.bcf.gz/_transcript.bcf.gz}
}
export -f select_transcript

function info_string {
    local filename=$1
    bcftools view -h ${filename} | grep '##INFO' | cut -f 3 -d = | cut -f 1 -d , | tr \\n , | sed 's/,$//'
}

function AF_filter {
    local AF_field=$1
    local threshold=$2
    echo "${AF_field}<=${threshold} || ${AF_field}>=$(echo 1-${threshold} | bc -l)"
}

function AF_filter_2 {
    local AF_field=$1
    local AF_field_2=$2
    local threshold=$3
    echo "($(AF_filter $AF_field $threshold)) || (${AF_field}=\".\" && ($(AF_filter $AF_field_2 $threshold)))"
}

function col_num {
    local filename=$1
    local column="$2"
    sed 1q $filename | tr \\t \\n | grep -n "${column}" | cut -f 1 -d :
}


###### Init ######
echo "Output folder: ${OUTPUT}, tmp folder: ${TMP}"
mkdir -p ${OUTPUT}
mkdir -p ${TMP}
cd ${OUTPUT}


###### Modules ######
module load samtools/1.11
module load parallel/20170222
module load bedtools/2.26.0
module load ensembl_api/89


###### Meta information ######
awk -F \\t '$4=="'${DOMAIN}'" {print $3}' ${INCLUDE} | sort -V > samples_${DOMAIN}.txt
${SCRIPT_DIR}/annotate_sample.py ${WD}/${CONF} > sample_annotation_${DOMAIN}.tsv


###### Transcripts coordinates in BED format ######
col_gene=$(col_num ${GENELIST} "${GL_gene}")
col_transcript=$(col_num ${GENELIST} "${GL_transcript}")
awk -F \\t -v OFS=\\t -v gene=${col_gene} -v transcript=${col_transcript} 'NR>1 {print $gene,$transcript}' ${GENELIST} > gene-transcript_${DOMAIN}.csv
perl ${SCRIPT_DIR}/get_genes_GRCh37.pl gene-transcript_${DOMAIN}.csv | sort -V > genes_${DOMAIN}.bed

# Add +/- 10bp flanks to the exons, merge and sort:
awk -F \\t -v OFS=\\t '{$2-=10; $3+=10; print $0}' genes_${DOMAIN}.bed | sort -V | bedtools merge -c 4,5 -o collapse | sort -V > genes_${DOMAIN}_sorted_merged.bed
# Get lists of transcripts per chromosome:
cut -f 1,4 genes_${DOMAIN}.bed | sort | uniq | awk -F \\t -v OFS=\\t '{print $2 > "transcripts_'${DOMAIN}'_chr"$1".txt"}'
CHROMOSOMES=$(cut -f 1 genes_${DOMAIN}.bed | sort -V | uniq)


###### Subset samples and regions ######
parallel --halt 1 -j ${CORES} --verbose "bcftools view -S samples_${DOMAIN}.txt -c1:nonmajor -R genes_${DOMAIN}_sorted_merged.bed -Ob -o ${TMP}/chr{}_${filename}_${DOMAIN}.bcf.gz ${MERGEDVCF}/chr{}_${filename}.bcf" ::: $CHROMOSOMES
# NB! -c1 vs -c1:nonmajor
parallel --halt 1 -j ${CORES} --verbose "bcftools index ${TMP}/chr{}_${filename}_${DOMAIN}.bcf.gz" ::: $CHROMOSOMES


###### Annotate with HGMD ######
remove_old=$(info_string ${TMP}/chr1_${filename}_${DOMAIN}.bcf.gz | tr , \\n | grep -E 'HGMD' | awk '{print "INFO/"$1}' | tr \\n , | sed s/,$//)
parallel --halt 1 -j ${CORES} --verbose "bcftools annotate -x ${remove_old} -Ob -o ${TMP}/chr{}_${filename}_${DOMAIN}_noHGMD.bcf.gz ${TMP}/chr{}_${filename}_${DOMAIN}.bcf.gz" ::: $CHROMOSOMES
parallel --halt 1 -j ${CORES} --verbose "bcftools index ${TMP}/chr{}_${filename}_${DOMAIN}_noHGMD.bcf.gz" ::: $CHROMOSOMES
HGMDfields="ID,$(info_string ${HGMD})"
echo "HGMD fields: $HGMDfields"
prefix_fields=$(prefix_anno HGMD ${HGMDfields})
parallel --halt 1 -j ${CORES} --verbose "bcftools annotate -c ${prefix_fields} -a ${HGMD} -Ob -o ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD.bcf.gz ${TMP}/chr{}_${filename}_${DOMAIN}_noHGMD.bcf.gz" ::: $CHROMOSOMES
rm ${TMP}/chr*_${filename}_${DOMAIN}_noHGMD.bcf.gz*
parallel --halt 1 -j ${CORES} --verbose "bcftools index ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD.bcf.gz" ::: $CHROMOSOMES


###### Annotate with ClinVar ######
CLNVfields="ID,$(info_string ${CLNV})"
echo "ClinVar fields: $CLNVfields"
prefix_fields=$(prefix_anno CLNV ${CLNVfields})
parallel --halt 1 -j ${CORES} --verbose "bcftools annotate -c ${prefix_fields} -a ${CLNV} -Ob -o ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD_CLNV.bcf.gz ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD.bcf.gz" ::: $CHROMOSOMES
parallel --halt 1 -j ${CORES} --verbose "bcftools index ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD_CLNV.bcf.gz" ::: $CHROMOSOMES


###### Split and subset transcripts from VEP-annotation ######
parallel --halt 1 -j ${CORES} --verbose "select_transcript {} ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD_CLNV.bcf.gz" ::: $CHROMOSOMES


###### Prioritise ######
# Can require customisation based on the project

# Step 1, AF and minOPR:
#  * keep gnomAD overall MAF <= 0.001 (same as AF<=0.001 OR AF>=0.999) OR, if in HGMD DM/DM? then MAF <= 0.025
#  * if missing in gnomAD, use AF_WGS10K with same logic
#  * remove minOPR<0.5 (other FILTERs are better than that)

filter_OPR="FILTER!=\"minOPR_0.5\""
filter_HGMD="HGMD_CLASS==\"DM,DM?\""
include="(${filter_OPR}) && \
(\
($(AF_filter_2 ${AF_field} ${AF_field_2} ${AF})) || \
(($(AF_filter_2 ${AF_field} ${AF_field_2} ${AF_HGMD})) && (${filter_HGMD}))\
)"
echo "Variant filtering 1, include: $include"

parallel --halt 1 -j ${CORES} --verbose "bcftools view -i '"$include"' -Ob -o ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD_CLNV_transcript_AF_minOPR.bcf.gz ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD_CLNV_transcript.bcf.gz" ::: $CHROMOSOMES
parallel --halt 1 -j ${CORES} --verbose "bcftools index ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD_CLNV_transcript_AF_minOPR.bcf.gz" ::: $CHROMOSOMES


# Step 2, consequence and impact:
#  * keep MODERATE|HIGH|splice_region_variant
#  * for minor ref (AF>0.5) or non-coding transcripts keep all non-synonymous
#  * on MT, remove upstream_gene_variant/downstream_gene_variant (because MT doesn't have any intronic/intergenic regions)

include1="VEP_IMPACT==\"MODERATE,HIGH\" || VEP_Consequence==\"splice_region_variant\""
include2="(${AF_field}>0.5 || (${AF_field}=\".\" && ${AF_field_2}>0.5) || VEP_BIOTYPE!=\"protein_coding\") && VEP_Consequence!=\"synonymous_variant\""
include="(${include1}) || (${include2})"
echo "Variant filtering 2, include: $include"

exclude="CHROM==\"MT\" && (VEP_Consequence==\"downstream_gene_variant,upstream_gene_variant\")"
echo "Variant filtering 2, exclude: $exclude"

parallel --halt 1 -j ${CORES} --verbose "bcftools view -i '"$include"' -Ou ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD_CLNV_transcript_AF_minOPR.bcf.gz | bcftools view -e '"$exclude"' -Ob -o ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD_CLNV_transcript_AF_minOPR_consequence.bcf.gz" ::: $CHROMOSOMES
parallel --halt 1 -j ${CORES} --verbose "bcftools index ${TMP}/chr{}_${filename}_${DOMAIN}_HGMD_CLNV_transcript_AF_minOPR_consequence.bcf.gz" ::: $CHROMOSOMES


# Step 3
# Sample level: remove allelic depth below 10

files_ordered=$(for i in {1..22} X Y MT; do fn=${TMP}/chr${i}_${filename}_${DOMAIN}_HGMD_CLNV_transcript_AF_minOPR_consequence.bcf.gz; if [[ -s $fn ]]; then echo $fn; fi; done)
echo "Merging results: ${files_ordered}"

format="$(echo [%${output_sample_format} ${output_variant_format}\\n] | sed 's/ /\\t%/g')"
echo "Format of the output table: $format"

# NB! | and & vs || and &&:
include1='((GT="alt" & ('$AF_field'<0.5 | ('$AF_field'="." | '$AF_field_2'<0.5))) | (GT!="AA" & GT!="A" & ('$AF_field'>0.5 | ('$AF_field'="." & '$AF_field_2'>0.5)))) & ((GT!="hap" & sSUM(FORMAT/AD)>=10) | (GT="hap" & sSUM(FORMAT/AD)>=5))'
# in case of homo-/hemizygous minor ref we won't have AD information. Include all of them:
include2='('$AF_field'>0.5 | ('$AF_field'="." & '$AF_field_2'>0.5)) & (GT=="RR" | GT=="R")'
include="($include1) | ($include2)"
echo "Sample filtering, include: $include"

echo ${output_sample_format} ${output_variant_format} | tr " " \\t > ${OUTPUT}/prioritised_variants_${DOMAIN}.tsv
# Not working without --naive:
bcftools concat --naive ${files_ordered} | bcftools query -f "$format" -i "$include" >> ${OUTPUT}/prioritised_variants_${DOMAIN}.tsv


# Tidy the result:
#  * split VEP_HGVSc and VEP_HGVSp
#  * in CLNV_CLNSIGCONF replace %3B with semicolon
#  * in HGMD_PHEN remove %2C (comma, probably an artifact)

${SCRIPT_DIR}/tidy.py ${OUTPUT}/prioritised_variants_${DOMAIN}.tsv > ${OUTPUT}/prioritised_variants_${DOMAIN}_tidy.tsv


###### Annotate (ethnicity, sex, consanguinity etc) ######
${SCRIPT_DIR}/merge.py -a ${OUTPUT}/prioritised_variants_${DOMAIN}_tidy.tsv \
    -b ${OUTPUT}/sample_annotation_${DOMAIN}.tsv -1 SAMPLE -2 SAMPLE -c SAMPLE -r 2 \
    > ${OUTPUT}/prioritised_variants_${DOMAIN}_tidy_samples.tsv


###### Annotate (genelist) ######
${SCRIPT_DIR}/merge.py -a ${OUTPUT}/prioritised_variants_${DOMAIN}_tidy_samples.tsv \
    -b ${GENELIST} -1 "VEP_Feature" -2 "${GL_transcript}" -c "VEP_EXON" \
    > ${OUTPUT}/prioritised_variants_${DOMAIN}_tidy_samples_genelist.tsv


exit 0
