# Contents: scripts used for manual investigation on INDEL (insertion and deletion) variants. 

    
# 1. Filtering for INDEL calls from total variant calls (SNP, INDEL, SV) from Strelka and Manta
    # bash script. bcftools. 
    # remove decoy sequences and filter for indels with more than one base pair. 
    
    #Strelka
    bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
     -f PASS \
     -i 'strlen(REF) - strlen(ALT) >= 2 || strlen(REF) - strlen(ALT) <= -2' \
        "Input_VCF_file" \
        -o indel_strelka.vcf
    
    
    #Manta __ small indels have CIGAR in format
    bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
     -f PASS -i 'INFO/CIGAR' \
      "Input_VCF_file" \
      -o indel_manta.vcf


# 2. Generating BED file from mutation calls (vcf) with flanking sequencing of 10bp
    # bash script. awk package.
  
      #!/bin/bash
      awk -v OFS='\t' '!/^#/{start=$2-10; end=$2+length($4)+10; if (start < 0) start = 0; print $1, start, end}' input_vcf_file > output_bed_file



# 3. Subsetting BAM file with BED file 
    # bash script. samtools package. input arguments defined. 
      #!/bin/bash
      
      # Check if the correct number of arguments is provided
      if [ "$#" -ne 3 ]; then
          echo "Usage: $0 <input_bed_file> <input_bam_file_name> <output_bam_file>"
          exit 1
      fi
      
      # Assign command-line arguments to variables
      input_bed_file="$1"
      input_bam_file_path="$2"
      output_bam_file="$3"
      
      # Path to samtools executable (adjust this based on your setup)
      samtools_path=""
      
      # Check if samtools executable exists and is executable
      if [ ! -x "$samtools_path" ]; then
          echo "Error: samtools not found or not executable"
          exit 1
      fi
      
      # Use samtools view to subset BAM file based on BED regions
      "$samtools_path" view -L "$input_bed_file" "$input_bam_file_path" -o "$output_bam_file"
      
      # Check if samtools view command was successful
      if [ "$?" -ne 0 ]; then
          echo "Error: samtools view command failed"
          exit 1
      fi
      
      echo "Subset BAM file created: $output_bam_file"

