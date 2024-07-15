# check integrity of all vcf files in the directory (arg1). 
    # used during debugging of internal mutation calling pipeline. 
    # output = .txt file of names of corrupted vcf files. 
    
#!/bin/bash

# Directory containing VCF files
VCF_DIR="path_to_your_directory"

# Output file for corrupted VCFs
OUTPUT_FILE="corrupted_vcfs.txt"

# Clear the output file if it already exists
> $OUTPUT_FILE

# Loop through all VCF files in the directory
for vcf_file in "$VCF_DIR"/*.vcf; do
    if [ -f "$vcf_file" ]; then
        # Check the VCF file for corruption
        bcftools view "$vcf_file" > /dev/null 2>&1
        if [ $? -ne 0 ]; then
            # If bcftools returns a non-zero exit code, the file is corrupted
            echo "$vcf_file" >> $OUTPUT_FILE
        fi
    fi
done

echo "Check completed. Corrupted VCF files (if any) are listed in $OUTPUT_FILE"
