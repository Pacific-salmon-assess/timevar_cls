#!/bin/bash -l
#
#SBATCH --array=0-863
#SBATCH --cpus-per-task=1
#SBATCH --job-name=samsim_rerun
#SBATCH --ntasks-per-node=1
#SBATCH --export=USER,LOGNAME,HOME,MAIL,PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#SBATCH --account=dfo_hpcmc_fm
#SBATCH --partition=standard
#SBATCH --time=28:00:00
#SBATCH --mem-per-cpu=6400M
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=catarinawor@gmail.com
#SBATCH --comment="image=registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu24.04"
#SBATCH --output=slurm_%a.out

cd /gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls/slurm_scripts  

export MAMBA_EXE=/home/spfm100/zhm001/miniforge3/bin/mamba
source /home/spfm100/zhm001/miniforge3/etc/profile.d/mamba.sh

mamba activate R_4.5.1

/gpfs/fs7/dfo/hpcmc/pfm/spfm100/zhm001/miniforge3/envs/R_4.5.1/bin/Rscript --vanilla parallel_run.R

mamba deactivate
