#!/bin/bash
#SBATCH -J pinetree
#SBATCH -n 2
#SBATCH -p normal
#SBATCH -o pinetree.o%j
#SBATCH -e pinetree.e%j
#SBATCH -t 18:00:00
#SBATCH -A A-bio7
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bjack913@gmail.com

module load launcher
module load gcc; module load python3

export LAUNCHER_WORKDIR=/work/03284/bjack/lonestar/phage_simulation/runs/012418


export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=/work/03284/bjack/lonestar/phage_simulation/runs/012418/commands_012418.sh

$LAUNCHER_DIR/paramrun
