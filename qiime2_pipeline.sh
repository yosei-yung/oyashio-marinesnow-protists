#!/bin/sh

OUT_DIR="/aptmp/ideas2/qingwei/KS-21/01_18S/"
MANIFEST="qiime_manifest_20220202_ks21.txt"
METADATA_FILE="KS-21_18S_metadata.tsv"

source /etc/profile.d/modules.sh

module load qiime2/2023.9
module load Python/3.9.5

# Step 1 -- import tools
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${OUT_DIR}${MANIFEST} \
  --input-format PairedEndFastqManifestPhred33 \
  --output-path ${OUT_DIR}out1-import.qza

Step2-quailty control with dada2
# use dada2 to trim off the primer sequence, to quailty control
qiime dada2 denoise-paired \
 --i-demultiplexed-seqs ${OUT_DIR}out1-import.qza \
 --p-trim-left-f 18 \
 --p-trim-left-r 20 \
 --p-trunc-len-f 260 \
 --p-trunc-len-r 240 \
 --o-denoising-stats ${OUT_DIR}out2_dada2_denoising-stats.qza \
 --o-representative-sequences ${OUT_DIR}out2-derep-sequences.qza \
 --o-table ${OUT_DIR}out2-derep-table.qza

# Step3-import reference data
# download lastest PR2 database from PR2 homepage

wget https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_mothur.fasta.gz
wget https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_mothur.tax.gz

# import sequence file
PR2_FILE="/aptmp/ideas2/qingwei/KS-21/01_18S/PR2v5/pr2_version_5.0.0_SSU_mothur.fasta"
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ${PR2_FILE} \
  --output-path ${OUT_DIR}out3-pr2_version_5.0.0_SSU_mothur.qza
# Import taxomony file
PR_FILE_TAX="/aptmp/ideas2/qingwei/KS-21/01_18S/PR2v5/pr2_version_5.0.0_SSU_mothur.tax"
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path ${PR_FILE_TAX} \
  --input-format HeaderlessTSVTaxonomyFormat \
  --output-path ${OUT_DIR}out3-pr2-ref-taxonomy.qza
# Step4-Training PR2 database used to classify 18S rRNA sequences
qiime feature-classifier extract-reads \
  --i-sequences ${OUT_DIR}out3-pr2_version_4.14.0_SSU_mothur.qza \
  --p-f-primer CYGCGGTAATTCCAGCTC \
  --p-r-primer AYGGTATCTRATCRTCTTYG \
  --p-min-length 200 \
  --p-max-length 500 \
  --o-reads ${OUT_DIR}out4-18Sref-seqsPR2.qza
# Naive Bayes Classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ${OUT_DIR}out4-18Sref-seqsPR2.qza \
  --i-reference-taxonomy ${OUT_DIR}out3-pr2-ref-taxonomy.qza \
  --o-classifier ${OUT_DIR}out4-PR2classifier.qza
# Step5-Classify 18S sequences
# Uisng q2-feature-classifier, containing several novel machine-learing and alignment-based for taxonomy classification
qiime feature-classifier classify-sklearn \
  --i-classifier ${OUT_DIR}out4-PR2classifier.qza \
  --i-reads ${OUT_DIR}out2-derep-sequences.qza \
  --o-classification ${OUT_DIR}out5-classify_PR2.qza

# Step6-filter singelton ASV/OTUs
qiime feature-table filter-features \
  --i-table ${OUT_DIR}out2-derep-table.qza \
  --p-min-frequency 2 \
  --m-metadata-file ${OUT_DIR}out5-classify_PR2.qza \
  --o-filtered-table ${OUT_DIR}out6-derep-clustered_ASV-final_table_no_singelton.qza
# Step7-Output table
#    a) ASV - table
qiime feature-table summarize \
  --i-table ${OUT_DIR}out6-derep-clustered_ASV-final_table_no_singelton.qza \
  --o-visualization ${OUT_DIR}out7-derep-clustered_ASV-final_table_no_singelton_OTU.qzv \
  --m-sample-metadata-file ${OUT_DIR}$METADATA_FILE
#    b) ASV - sequence list
qiime feature-table tabulate-seqs \
  --i-data ${OUT_DIR}out2-derep-sequences.qza \
  --o-visualization ${OUT_DIR}out7.2-seqs-nonchimeras_clustered_ASV_OTU_Seqs.qzv
#    c) taxonomic read table as tsv and   ##   f) OTU table with reads as tsv
mkdir -p ${OUT_DIR}export_for_ASV_table
qiime tools export \
 --input-path ${OUT_DIR}out7-derep-clustered_ASV-final_table_no_singelton_OTU.qzv  \
 --output-path ${OUT_DIR}export_for_ASV_table

/bin/mv -f ${OUT_DIR}export_for_ASV_table/feature-table.biom \
${OUT_DIR}export_for_ASV_table/out7-derep-clustered_ASV-final_table_no_singelton_feature-table.biom
biom convert \
  -i ${OUT_DIR}export_for_ASV_table/out7-derep-clustered_ASV-final_table_no_singelton_feature-table.biom \
  -o ${OUT_DIR}export_for_ASV_table/out7-derep-clustered_ASV-final_table_no_singelton_OTUs.tsv \
  --to-tsv
#   d) OTU table with assigned taxonomy
# exporting classification ASV-tax list from qiime
qiime tools export \
 --input-path ${OUT_DIR}out5-classify_PR2.qza \
 --output-path ${OUT_DIR}export_tax_ASV_list
# exporting ASV table without annotation from qiime
qiime tools export \
 --input-path ${OUT_DIR}out2-derep-table.qza\
 --output-path ${OUT_DIR}export_tax_ASV_list
# Step8- Rename the classification (ASV - tax) list header"
# replace fist line of taxonomy.tsv with "#ASVID       taxonomy        confidence"
awk 'BEGIN {FS="\t"; OFS="\t"; print "#OTUID", "taxonomy", "confidence"} NR > 1' ${OUT_DIR}export_tax_ASV_list/taxonomy.tsv ${OUT_DIR}export_tax_ASV_list/taxonomy.tsvASV_renamed.tsv
# a) adding metadata (taxonomy) to the ASV table (feature-table.biom)
biom add-metadata \
  -i ${OUT_DIR}export_tax_ASV_list/feature-table.biom \
  --observation-metadata-fp ${OUT_DIR}export_tax_ASV_list/taxonomy.tsvASV_renamed.tsv \
  --sc-separated taxonomy \
  -o ${OUT_DIR}export_tax_ASV_list/table-with-taxonomy_ASV.biom
# b)exporting the ASV table (feature-table.biom) WITH the added classification (tax)
biom convert \
 -i ${OUT_DIR}export_tax_ASV_list/table-with-taxonomy_ASV.biom \
 --to-tsv \
 --header-key taxonomy \
 -o ${OUT_DIR}export_tax_ASV_list/out8-derep-filtered_table_clustered_ASV_with_tax.tsv


