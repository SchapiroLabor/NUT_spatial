#!/bin/bash
#SBATCH --job-name=cellcharter
#SBATCH --partition=cpu-single
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=300gb
#SBATCH --array=4

module load devel/miniconda/3
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate cellcharter_env

python cellcharter_batch.py
