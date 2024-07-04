# script used to convert BAM file back to FASTQ file (paired-read)
#  sort both bam and fastq files.

# required argument: sample ID, path to sample file 



# Ensure that an input sample is provided

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_sample>"
    exit 1
fi

# Get the input sample
input_sample=$1

# Define the directory containing the BAM file
bam_dir="$2/${input_sample}"

# Ensure the directory exists
if [ ! -d "$bam_dir" ]; then
    echo "Directory $bam_dir does not exist!"
    exit 1
fi

# Find the BAM file in the directory
bam_files=("$bam_dir"/*.bam)

# Check if exactly one BAM file is found
if [ "${#bam_files[@]}" -ne 1 ]; then
    echo "Expected one BAM file in $bam_dir, but found: ${#bam_files[@]}"
    exit 1
fi

# Get the input BAM file path
input_bam="${bam_files[0]}"

# Copy the BAM file to the current launch directory & sort bam 
cp "$input_bam" .
samtools sort "$input_bam"

# Get the BAM file name
bam_file_name=$(basename "$input_bam")

# Extract prefix from BAM file name (without extension)
prefix=${bam_file_name%%.bam}

# Convert BAM to FASTQ using samtools
samtools fastq -@ ${SLURM_CPUS_PER_TASK} -1 "${prefix}_R1.fastq" -2 "${prefix}_R2.fastq" "$bam_file_name"

# Check if samtools command was successful
if [ $? -ne 0 ]; then
    echo "samtools fastq failed!"
    exit 1
fi

# Function to sort and compress FASTQ files
sort_and_compress_fastq() {
    local fastq_file=$1
    local sorted_fastq_file="sorted_${fastq_file}.gz"

    awk '{ 
        if (NR % 4 == 1) { 
            header = $0; 
            getline seq; 
            getline qheader; 
            getline qual; 
            print header "\t" seq "\t" qheader "\t" qual 
        } 
    }' "$fastq_file" | sort -k1,1 | awk -F"\t" '{ 
        print $1 "\n" $2 "\n" $3 "\n" $4 
    }' | gzip > "$sorted_fastq_file"

    # Check if the sorting and compression were successful
    if [ $? -ne 0 ]; then
        echo "Failed to sort and compress $fastq_file"
        exit 1
    fi

    # Remove the intermediate FASTQ file
    rm -f "$fastq_file"
    if [ $? -ne 0 ]; then
        echo "Failed to remove intermediate file $fastq_file"
        exit 1
    fi
}

# Sort and compress R1 FASTQ file
sort_and_compress_fastq "${prefix}_R1.fastq"

# Sort and compress R2 FASTQ file
sort_and_compress_fastq "${prefix}_R2.fastq"

# Remove the copied BAM file
rm -f "$bam_file_name"
if [ $? -ne 0 ]; then
    echo "Failed to remove copied BAM file $bam_file_name"
    exit 1
fi

echo "Processing complete. Sorted and compressed files:"
echo "sorted_${prefix}_R1.fastq.gz"
echo "sorted_${prefix}_R2.fastq.gz"
