#!/bin/bash            
                                                        
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=3882
#SBATCH --time=5:59:00

TMPDIR=/data/scratch/${USER}/mcrcache_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR
export MCR_CACHE_ROOT=$TMPDIR
echo "$SLURM_ARRAY_TASK_ID"
srun ./run_OOS_slurm_IF.sh /opt/MATLAB/R2017a $SLURM_ARRAY_TASK_ID
rm -rf $TMPDIR