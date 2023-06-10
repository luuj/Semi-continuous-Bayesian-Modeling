#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 6-00:00
#SBATCH --mem 12G

module load JAGS/4.3.0-fasrc01
module load R/4.1.0-fasrc01
export R_LIBS_USER=$HOME/apps/4.1.0:$R_LIBS_USER
Rscript main.R ${PARAM1}