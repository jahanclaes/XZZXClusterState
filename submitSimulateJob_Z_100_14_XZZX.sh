#!/bin/bash
#SBATCH --job-name=Cluster_Z_100_14_XZZX
#SBATCH --requeue
#SBATCH --partition scavenge
#SBATCH --time=24:00:00

#SBATCH --mem=20G

module load matplotlib/3.3.4-foss-2020b
for i in {1..1000}; do python SimulateThresholds.py 100 0.062 14 ; done