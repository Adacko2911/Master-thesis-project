# Load modules
module load 2021
module load R/4.1.0-foss-2021a
module load Stopos/0.93-GCC-10.3.0

# Make the script executable
chmod +x $HOME/convergence_run/mhmm/doit

# Submit jobs
sbatch $HOME/convergence_run/mhmm/jobs/job1

