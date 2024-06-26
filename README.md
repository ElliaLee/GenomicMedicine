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
