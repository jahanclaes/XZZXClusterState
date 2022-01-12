#!/bin/bash
#SBATCH --job-name=Cluster_Z_5000_7_XZZX
#SBATCH --requeue
#SBATCH --partition scavenge
#SBATCH --time=24:00:00

module load matplotlib/3.3.4-foss-2020b
for i in {1..1000}; do python SimulateThresholds.py 5000 0.081 7 ; done