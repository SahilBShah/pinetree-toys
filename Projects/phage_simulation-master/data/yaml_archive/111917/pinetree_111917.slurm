#!/bin/bash
#SBATCH -J pinetree_ko
#SBATCH -n 3
#SBATCH -p normal
#SBATCH -o pinetree_ko.o%j
#SBATCH -e pinetree_ko.e%j
#SBATCH -t 12:00:00
#SBATCH -A A-bio7
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bjack913@gmail.com

module load launcher
module load gcc; module load python3

export LAUNCHER_WORKDIR=/work/03284/bjack/lonestar/phage_simulation/runs/111917


export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=/work/03284/bjack/lonestar/phage_simulation/runs/111917/commands_111917.sh

$LAUNCHER_DIR/paramrun
