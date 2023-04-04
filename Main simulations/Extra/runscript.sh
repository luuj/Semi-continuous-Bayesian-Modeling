#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 2-00:00
#SBATCH --mem 8G
#SBATCH --mail-type=END
#SBATCH --mail-user=jluu@g.harvard.edu

module load R/4.1.0-fasrc01
export R_LIBS_USER=$HOME/apps/4.1.0:$R_LIBS_USER
Rscript main.R $1