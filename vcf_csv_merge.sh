# written to create mutation list in csv format from multiple vcf files to be used for INVAR2 pipeline. 
# argument needed = directory path of vcf files.

#!/bin/bash

# Directory containing the VCF files
vcf_dir=""
bcftools=""

# Directory for output CSV files
output_dir="./"
mkdir -p "$output_dir"

# Process each VCF file
for vcf_file in "$vcf_dir"/*.vcf; do
  # Extract the base name of the file without extension
  base_name=$(basename "$vcf_file" .vcf)
  
  # Extract patientID from the filename using sed
  patient_id=$(echo "$base_name" | sed -n 's/^filtered_\([^\.]*\)\..*/\1/p')

  # Output file path
  output_file="$output_dir/$patient_id.csv"
  
  # Extract information
  ${bcftools} query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF\t]\n' "$vcf_file" > "$output_dir/raw.txt"

  # Only get tumour AF and add patient information
  awk -v patient_id="$patient_id" -F'\t' 'BEGIN {OFS=","} {
    print $1, $2, $3, $4, $5, patient_id
  }' "$output_dir/raw.txt" > "$output_file"

  # Remove the temporary raw file
  rm "$output_dir/raw.txt"
done

# Create a merged CSV file
merged_file="$output_dir/MutationCalls.csv"

# Add the header to the merged file
echo 'CHROM,POS,REF,ALT,TUMOUR_AF,PATIENT' > "$merged_file"

# Append each CSV file to the merged file
for csv_file in "$output_dir"/*.csv; do
  cat "$csv_file" >> "$merged_file"
done

echo "CSV files merged into: $merged_file"
