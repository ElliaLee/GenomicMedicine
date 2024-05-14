"""
  Python scripts written for debugging INVAR pipeline

  Contents:
      1. Subsetting mutation calls of CSV file
      2. Converting mutation list in txt file format into CSV file with correct format
      3. Filtering sample list file that corresponse to available samples
  """

# 1. Subsetting mutation calls (csv file) to test if error is due to memory issue
      # for debugging error status 140. 

    import os
    import pandas as pd
    
    directory = '${OutputDirectory}'
    os.chdir(directory)
    
    
    mutations_df = pd.read_csv('${MutationCallFile}')
    mut_onChr1 = mutations_df[mutations_df['CHROM'] == 'chr1']
    
    mut_onChr1.to_csv('MutationCalls_CHR1', index=False)

# 2. Converting filtered mutation calls in .txt file to .csv file with correct format
      # for debugging errors that indicate "patients do not have mutations" that arise due to incorrect format of mutation list

    import pandas as pd
    import os
    
    folderPath = "${DirectoryToMutationListFile}"
    filePath = folderPath
    
    fileList = [file for file in os.listdir(filePath) if file.startswith("${StudyName}")]
    
    mainFile_list = []
    
    for file_name in fileList:
        file_path = os.path.join(filePath, file_name)
        file = pd.read_csv(file_path, delimiter="\t")
        
        # Overwrite the patient ID number from patient_XXX to just XXX
        file["PATIENT"] = file_name[3:7]
        
        mainFile_list.append(file)
    
    mainFile = pd.concat(mainFile_list, ignore_index=True)

    # Renaming column names to be correctly recognised by INVAR
    mainFile = mainFile.rename(columns={"Chromosome":"CHROM", "Position" : "POS" , "Reference allele" : "REF", "Alternate allele":"ALT", "Allele fraction (tumour)":"TUMOUR_AF"})
    
    mainFile.to_csv(os.path.join(folderPath, "MutationCalls.csv"), index=False)


# 3. Subsetting total sample list file (csv) that only includes available sample. 
    # Used during pre-processing for alignment, when sample sheet includes samples that are not sequenced. 
      # txt file of samples that have been sequenced is created prior
      import csv
      import os
      
      def filter_nonexistent_files(csv_file_path, existing_files_path):
          with open(existing_files_path, 'r') as existing_file:
              existing_files = {line.strip() for line in existing_file if line.strip()}
          
          filtered_rows = []
          with open(csv_file_path, 'r') as csvfile:
              csv_reader = csv.reader(csvfile)
              header = next(csv_reader)
              filtered_rows.append(header)
      
              for row in csv_reader:
                  if len(row) >1:
                      filename_col0 = row[0].strip()
                      filename_col1 = row[1].strip()
                      if filename_col0 in existing_files or filename_col1 in existing_files:
                          filtered_rows.append(row)
      
          # Rewrite the CSV file with only the filtered rows
          with open(csv_file_path, 'w', newline='') as csvfile:
              csv_writer = csv.writer(csvfile)
              csv_writer.writerows(filtered_rows)
      
      filter_nonexistent_files('${pathToSampleSheet_csv}', '${PathToListofAvailableSample_txt}')

