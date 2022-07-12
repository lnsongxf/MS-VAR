#!/bin/bash            
                                                        
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=3882
#SBATCH --time=5:59:00

# path to MATLAB runtime
MCR=/apps/MATLAB/R2017a
 
# create an independent MATLAB runtime cache for each instance,
# since they will otherwise try to use the same path
MCR_CACHE_ROOT=/lcl/data/lscratch/${USER}/mcrcache_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export MCR_CACHE_ROOT

# create cache dir
mkdir -p "$MCR_CACHE_ROOT"
echo "$SLURM_ARRAY_TASK_ID"

cd `pathf $SLURM_SUBMIT_DIR`
srun run_OOS_slurm_RS.sh $MCR $SLURM_ARRAY_TASK_ID

# delete cache dir
rm -rf "$MCR_CACHE_ROOT"