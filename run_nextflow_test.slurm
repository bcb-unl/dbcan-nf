#!/bin/bash
#SBATCH --job-name=nf-dbcan-test
#SBATCH --output=logs/nf-test-%j.out
#SBATCH --error=logs/nf-test-%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=batch
#SBATCH --qos=ac_yourgroup   # Uncomment and replace with your group name if you have acknowledgement credit

# Create logs directory if it doesn't exist
mkdir -p logs

# Load necessary modules
# Adjust the following commands according to your cluster's module system
module load nextflow
module load singularity   # If using Singularity/Apptainer

# Set environment variables
# Nextflow work directory - use work storage instead of home
export NXF_WORK=/work/$USER/nf-work
mkdir -p $NXF_WORK

# Singularity/Apptainer cache directory
export NXF_SINGULARITY_CACHEDIR=/work/$USER/.singularity_cache
export SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR
mkdir -p $NXF_SINGULARITY_CACHEDIR

# Temporary directory - use scratch space if available
export TMPDIR=/scratch/$USER/nf_test_tmp
mkdir -p $TMPDIR

# Set Java options (required by Nextflow)
export NXF_OPTS="-Xms1g -Xmx4g"

# Change to project directory
cd /array1/xinpeng/dbcan_new_dev/dbcan-nf

# Run Nextflow pipeline test
# Using test profile and singularity profile
nextflow run main.nf \
    -profile test,singularity,slurm \
    --input dbcan-nf-test_data/testsamplesheet.csv \
    --outdir ./test_results \
    -c conf/slurm.config \
    -resume

# Check exit status
if [ $? -eq 0 ]; then
    echo "Nextflow pipeline completed successfully!"
else
    echo "Nextflow pipeline failed. Check error logs."
    exit 1
fi
