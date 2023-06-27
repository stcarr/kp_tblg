#!/bin/bash
#SBATCH -J KOTBH
#SBATCH -p kaxiras,shared
#SBATCH -t 160:00:00
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -o KP_testBZ.out
#SBATCH -e KP_testBZ.err
#SBATCH --mem=100000
 
module load matlab/R2017a-fasrc02
srun -n 1 -c 2 matlab-default -nosplash -nodesktop -r "global_vars; run('/home/stc/devspace/codes/shiang_kp_relaxed_tblg/kp_construction_tools/parameter_extract_job.m'); exit"


