#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name="Odin5"
#SBATCH -o errors/%x.%j.%a.out           # output
#SBATCH -e errors/%x.%j.%a.err           # errors
#SBATCH --partition=Lake
#SBATCH --mem=1GB
#SBATCH --time=6-00:00:00
#SBATCH --array=0-59
# with:
# %N = nodename
# %x = jobname
# %j = jobid
# %a = arrayid

jump_num=$((($SLURM_ARRAY_TASK_ID%10)+1))
rep_num=$(($SLURM_ARRAY_TASK_ID/10))

echo $jump_num
echo $rep_num

module load R

Rscript --vanilla run_odin_m5_twomethods_jumps_reps_lrompsi.R Influenza RSV newvalidateSimpleInter_ append cyc $jump_num all $rep_num 5000000 50
