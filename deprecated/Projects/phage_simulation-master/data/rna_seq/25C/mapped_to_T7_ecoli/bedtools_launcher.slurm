#!/bin/bash
#SBATCH -J multicov
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p normal
#SBATCH -o multicov.o%j
#SBATCH -e multicov.e%j
#SBATCH -t 02:00:00
#SBATCH -A A-bio7

module load launcher
module load samtools
module load bedtools




export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=count_commands.sh

$LAUNCHER_DIR/paramrun
