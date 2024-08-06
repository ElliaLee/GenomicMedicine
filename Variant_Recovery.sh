# Written to compare mutation calling of tumour_vs_buffyCoat and tumour_vs_plasma 
# Interested to find out if actual mutations were called using plasma, but were filtered out 

#!/bin/bash
BASE_ID=$1
tum_cfDNA_path=""
tum_buffy_path=""
bcftools=""
tools=""

mkdir ${BASE_ID}
cd ${BASE_ID}

# get tumour vs plasma mutation calls (vcf files: .refiltered.pass.clean.vcf)
cp ../${tum_cfDNA_path}/${BASE_ID}.snv.metrics.re_filtered.pass.clean.vcf ./
mv ${BASE_ID}.snv.metrics.re_filtered.pass.clean.vcf ${BASE_ID}.plasma.mutcalled.vcf

# get tumour vs buffy mutation calls (vcf file)
cp ../${tum_buffy_path}/${BASE_ID}.snv.metrics.re_filtered.pass.clean.vcf ./
mv ${BASE_ID}.snv.metrics.re_filtered.pass.clean.vcf ${BASE_ID}.buffy.mutcalled.vcf

# clean (remove decoy sequences), gzip, index input vcf files
grep -E '^#|^chr[1-9]\b|^chr1[0-9]\b|^chr2[0-2]\b|^chrX\b|^chrY\b' ${BASE_ID}.plasma.mutcalled.vcf > ${BASE_ID}.plasma.vcf
${tools}/bgzip ${BASE_ID}.plasma.vcf
${tools}/tabix -p vcf ${BASE_ID}.plasma.vcf.gz

grep -E '^#|^chr[1-9]\b|^chr1[0-9]\b|^chr2[0-2]\b|^chrX\b|^chrY\b' ${BASE_ID}.buffy.mutcalled.vcf > ${BASE_ID}.buffy.vcf
${tools}/bgzip ${BASE_ID}.buffy.vcf
${tools}/tabix -p vcf ${BASE_ID}.buffy.vcf.gz

# get overlapping mutations called by both plasma & buffy / output VCF will have info from plasma vcf 
${bcftools} isec -n=2 -w 1 ${BASE_ID}.plasma.vcf.gz ${BASE_ID}.buffy.vcf.gz -O z -o ${BASE_ID}_overlap.vcf.gz
${tools}/tabix -p vcf ${BASE_ID}_overlap.vcf.gz
overlap_count=$(${bcftools} view --no-header ${BASE_ID}_overlap.vcf.gz | wc -l)

# get plasma specific mutation calls 
${bcftools} isec -C ${BASE_ID}.plasma.vcf.gz ${BASE_ID}_overlap.vcf.gz -w 1 -O z -o ${BASE_ID}_plasma_specific.vcf.gz
${tools}/tabix -p vcf ${BASE_ID}_plasma_specific.vcf.gz
plasma_specific_count=$(${bcftools} view --no-header ${BASE_ID}_plasma_specific.vcf.gz | wc -l)

# get buffy specific mutation calls
${bcftools} isec -C ${BASE_ID}.buffy.vcf.gz ${BASE_ID}_overlap.vcf.gz -w 1 -O z -o ${BASE_ID}_buffy_specific.vcf.gz
${tools}/tabix -p vcf ${BASE_ID}_buffy_specific.vcf.gz
buffy_specific_count=$(${bcftools} view --no-header ${BASE_ID}_buffy_specific.vcf.gz | wc -l)

# get overlapping (match)% = overlap / (plasma specific + buffy specific + overlap ) * 100
all_count=$(awk "BEGIN {print($overlap_count + $plasma_specific_count + $buffy_specific_count)}")
match_percentage=$(awk "BEGIN {print(${overlap_count}/${all_count} * 100)}")
echo "Match %: $match_percentage % match between mutation calls from tumour_vs_buffy and tumour_vs_plasma" > results.txt
echo "Overlap count = $overlap_count" >> results.txt
echo "Buffy specific count = $buffy_specific_count . vcf file saved as ${BASE_ID}_buffy_specific.vcf.gz" >> results.txt
echo "Plasma specific count = $plasma_specific_count . vcf file saved as ${BASE_ID}_plasma_specific.vcf.gz" >> results.txt

# try to find missing calls (called by buffy but not by plasma) in pre-filtered vcf files: re_filtered.clean.vcf & .vcf 
cp ../${tum_cfDNA_path}/${BASE_ID}.snv.metrics.re_filtered.clean.vcf ./
mv ${BASE_ID}.snv.metrics.re_filtered.clean.vcf ${BASE_ID}.plasma.re_filtered.clean.mutlist.vcf

cp ../${tum_cfDNA_path}/${BASE_ID}.vcf ./
mv ${BASE_ID}.vcf ${BASE_ID}.plasma.pre_filter.mutlist.vcf


# clean (remove decoy sequences), gzip, index  vcf files
grep -E '^#|^chr[1-9]\b|^chr1[0-9]\b|^chr2[0-2]\b|^chrX\b|^chrY\b' ${BASE_ID}.plasma.re_filtered.clean.mutlist.vcf > ${BASE_ID}.plasma.re_filtered.clean.vcf
${tools}/bgzip ${BASE_ID}.plasma.re_filtered.clean.vcf
${tools}/tabix -p vcf ${BASE_ID}.plasma.re_filtered.clean.vcf.gz

grep -E '^#|^chr[1-9]\b|^chr1[0-9]\b|^chr2[0-2]\b|^chrX\b|^chrY\b' ${BASE_ID}.plasma.pre_filter.mutlist.vcf > ${BASE_ID}.plasma.pre_filter.vcf
${tools}/bgzip ${BASE_ID}.plasma.pre_filter.vcf
${tools}/tabix -p vcf ${BASE_ID}.plasma.pre_filter.vcf.gz


# find mutation calls lost during filtering (not in final plasma vcf file, but called by buffy) / output include info from plasma vcf file 
${bcftools} isec -n=2 -w 1 ${BASE_ID}.plasma.re_filtered.clean.vcf.gz ${BASE_ID}_buffy_specific.vcf.gz -O z -o ${BASE_ID}_overlap.re_filtered.clean.vcf.gz
${tools}/tabix -p vcf ${BASE_ID}_overlap.re_filtered.clean.vcf.gz
overlap_re_filtered_count=$(${bcftools} view --no-header ${BASE_ID}_overlap.re_filtered.clean.vcf.gz | wc -l)

${bcftools} isec -n=2 -w 1 ${BASE_ID}.plasma.pre_filter.vcf.gz ${BASE_ID}_buffy_specific.vcf.gz -O z -o ${BASE_ID}_overlap.pre_filter.vcf.gz
${tools}/tabix -p vcf ${BASE_ID}_overlap.pre_filter.vcf.gz
overlap_pre_filter_count=$(${bcftools} view --no-header ${BASE_ID}_overlap.pre_filter.vcf.gz | wc -l)

# from pre-filter.vcf calls, get mutations filtered ONLY from default mutect2 filter (exclude re_filtered.clean.vcf calls)
${bcftools} isec -C -w 1 ${BASE_ID}_overlap.pre_filter.vcf.gz ${BASE_ID}_overlap.re_filtered.clean.vcf.gz -O z -o ${BASE_ID}_overlap.only_by_default.vcf.gz
${tools}/tabix -p vcf ${BASE_ID}_overlap.only_by_default.vcf.gz
filtered_only_by_default=$(${bcftools} view --no-header ${BASE_ID}_overlap.only_by_default.vcf.gz | wc -l)


# mutations lost during snv.metrics (no pass)
echo " " >> results.txt
echo "total # of mutations lost during filtering = $overlap_pre_filter_count. vcf file saved as ${BASE_ID}_overlap.pre_filter.vcf.gz" >> results.txt
echo "filtered out mutations by default mutect2 filtering = $filtered_only_by_default . vcf file saved as ${BASE_ID}_overlap.only_by_default.vcf.gz" >> results.txt
echo "filtered out mutations by snv.metrics (no pass) = $overlap_re_filtered_count . vcf file saved as ${BASE_ID}_overlap.re_filtered.clean.vcf.gz" >> results.txt

# proportion of normal_artefacts 
echo " " >> results.txt
echo "proportion of actual mutations flagged as normal_artefact:" >> results.txt

pre_filter_na_count=$(${bcftools} view --no-header -f "normal_artifact" ${BASE_ID}_overlap.pre_filter.vcf.gz | wc -l )
pre_filter_na_percentage=$(awk "BEGIN {print(${pre_filter_na_count}/${overlap_pre_filter_count} * 100)}")
echo "normal_artifact % pre-filtering =$pre_filter_na_percentage %" >> results.txt


re_filter_na_count=$(${bcftools} view --no-header -f "normal_artifact" ${BASE_ID}_overlap.re_filtered.clean.vcf.gz | wc -l )
re_filter_na_percentage=$(awk "BEGIN {print(${re_filter_na_count}/${overlap_re_filtered_count} * 100)}")
echo "normal_artifact % after snv_metrics filter(no pass) = $re_filter_na_percentage %" >> results.txt

pass_na_count=$(${bcftools} view --no-header -f "normal_artifact" ${BASE_ID}.plasma.vcf.gz | wc -l )
echo "after all the filtering, plasma mutList had $pass_na_count normal_artefacts" >> results.txt
