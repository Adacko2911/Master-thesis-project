# Load modules
module load 2021
module load R/4.1.0-foss-2021a
module load Stopos/0.93-GCC-10.3.0

# Check pools
stopos pools

# Delete previous poolsets
stopos purge -p poolset1

# Create pool
stopos create -p poolset1 
stopos add $HOME/mhmm_run/parameters/parm_set1 -p poolset1

# Make the script executable
chmod +x $HOME/mhmm_run/doit

# Submit jobs
sbatch $HOME/mhmm_run/jobs/job1


