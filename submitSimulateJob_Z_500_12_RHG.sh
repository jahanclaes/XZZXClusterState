#!/bin/bash
#SBATCH --job-name=Cluster_Z_500_12_RHG
#SBATCH --requeue
#SBATCH --partition scavenge
#SBATCH --time=24:00:00

module load matplotlib/3.3.4-foss-2020b
for i in {1..1000}; do python SimulateThresholds.py 500 0.03 12 RHG ; done