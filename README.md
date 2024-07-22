This repository exhibits projects from Genomic Medicine Master's Course at University of Cambridge.

Contents:

1. ClinicalData_statistics_MachineLearning.sh
    - R-script
    - Supervised by Dr.Hearn (Academic Tutor)
    - Project involved:
        - designing data structure
        - processing clinical data
        - performing statistical test
        - building machine learning models

     
2. Indel_Calling.sh
    - Bash scripts
    - Part of Research Project based in CRUK Cambridge Institute.
    - Used for manual INDEL Investigation.
    - includes scripts for:
        - Caling INDEl calls on VCF files from Strelka and Manta    
        - Generating BED file with final VCF file 
        - Subsetting BAM file with INDEL BED file to facilitate visualisation of IGV. 

3. INVAR_debug.py
    - Python scripts
    - Part of Research Project based in CRUK Cambridge Institute
    - Used for debugging errors encountered during running INVAR pipeline
    - includes sciprts for:
        - Subsetting mutation lists (for testing peak memory usage)
        - Generating mutation list (csv file) in correct format from .txt files obtained from multiple patients
        - Subsetting total sample sheet for only those with available sequence data

4. bamToFastq.sh
   - Bash script
   - Part of Research Project based in CRUK Cambridge Institute
   - Used to convert BAM file back to FASTQ
   - written to realign reads using internal pipeline.
  
5. check_vcf.sh
    - bash script
    - part of Rearch Project based in CRUK Cambridge Institute.
    - Written to check integrity of all vcf files in the directory & save names of corrupted files in a txt file. 
    - Used to debug internal mutation calling pipeline. 

6. vcf_csv_merge.sh
   - bash script
   - part of research project based in CRUK Cambridge Institute.
   - written to create a merged mutation list (.csv) file from multiple vcf files, for INVAR2 analysis.
  
7. vcf_txt.sh
   - bash script
   - part of research project based in CRUK Cambridge Institute.
   - written to debug merged vcf files having duplicate of same samples.
   - output files include: mutation_list.txt and allele_count.txt
   - mutation list.txt can be used for INVAR2 analysis.
   - allele count.txt can be used to check allele count of the variant in normal sample (filtered for allele count >1)
   
