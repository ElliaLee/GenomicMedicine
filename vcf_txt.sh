# written to overcome merged vcf files having duplicate of sample info. 

# arguments needed: patient ID, merged vcf file name 
    # merged vcf file has a duplicate of tumor & normal samples

# aim: create patientID.txt file with following extracted info:
    # CHROM / POS / REF / ALT  / TUMOUR AF
    # depending on which AF was chosen during merging of AF, one tumour sample will be 'empty' ('.')
    # extract tumour AF, from the duplicated tumour sample 

# additionally: create "${patientID}_allele_count.txt" to check allele count of the variant in normal sample and filter for those with at least one allele count
  # selectively extract "INFO/AD" from both normal samples
  # extract allele count (number after 'AD,')
  # filter for variants with allele count > 1 



bcftools=""

# arguments
patientID=$1
cd ./${patientID}
vcf_file=""

# get sample names 
${bcftools} query -l ${vcf_file} > sample_names.txt                                        
tum_1=$(head -n 2 sample_names.txt | tail -n 1)  
tum_2=$(tail -n 1 sample_names.txt)              

# extract to txt file
${bcftools} query -s ${tum_1},${tum_2} -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF\t]\n' ${vcf_file} > raw.txt 


# merge tumour AF 
awk -F'\t' '{
    if ($5 ~ /^[0-9]+\.[0-9]+$/) {
        $7 = $5
    } else if ($6 ~ /^[0-9]+\.[0-9]+$/) {
        $7 = $6
    }
    $5 = $6 = ""
    OFS=","; print $1, $2, $3, $4, $7
}' raw.txt > raw2.txt

awk -F',' -v new_col="${patientID}" 'BEGIN {OFS=FS} {print $0, new_col}' raw2.txt > "${patientID}_mutations.csv"

sed -i '1iCHROM,POS,REF,ALT,TUMOUR_AF,PATIENT' "${patientID}_mutations.csv"

# delete intermediate file
rm raw*




# for variants with normal allele count >= 1
normal_1=$(head -n 1 sample_names.txt) 
normal_2=$(tail -n 2 sample_names.txt | head -n 1)

# extract info
${bcftools} query -s ${normal_1},${normal_2} -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\t[%AF\t]\n' ${vcf_file} > raw.txt

awk -F'\t' '{
    if ($5 == ".") {
        $10 = $6
    } else if ($6 == "."){
        $10 = $5
    }
    OFS="\t"; print $1, $2, $3, $4, $8, $9, $10
}' raw.txt > raw2.txt

awk -F'\t' '{
    if ($5 ~ /^[0-9]+\.[0-9]+$/) {
        $8 = $5
    } else if ($6 ~ /^[0-9]+\.[0-9]+$/) {
        $8 = $6
    }
    $5 = $6 = ""
    OFS="\t"; print $1, $2, $3, $4, $7, $8
}' raw2.txt > raw3.txt

awk -F'\t' '{
    split($5, values, ",")

    # Extract the second value
    second_value = values[2]

    # Only print the row if the second value is not 0
    if (second_value != "0") {
        $5 = second_value
        OFS="\t"
        print $1, $2, $3, $4, $6, $5
    }
}' raw3.txt > "${patientID}_allele_count.txt"

sed -i '1iCHROM\tPOS\tREF\tALT\tNORMAL_AF\tALLELE_COUNT' "${patientID}_allele_count.txt"


rm raw*
rm sample_names.txt

mkdir mutation_list
mv "${patientID}_mutations.csv" ./mutation_list
mv "${patientID}_allele_count.txt" ./mutation_list
    
