# Load modules
module load 2021
module load R/4.1.0-foss-2021a
module load Stopos/0.93-GCC-10.3.0

# Check pools
stopos pools

# Delete previous poolsets
stopos purge -p poolset5
# Create pool
stopos create -p poolset5
stopos add $HOME/medhmm_run/parameters/parm_set5 -p poolset5

# Make the script executable
chmod +x $HOME/medhmm_run/doit5

# Submit jobs
sbatch $HOME/medhmm_run/jobs/job5


