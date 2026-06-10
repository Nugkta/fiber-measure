#!/bin/bash
# Full-dataset fibre measurement run.
#
# Dependencies: uv project at /net/scratch/j56806hx/spins-cv/fibrecv (fibrecv package).
# Inputs:  all "masp2 *_*.jpg" under /net/scratch/j56806hx/spins-cv/Images MasP2.
# Output:  per-image artifacts under /net/scratch/j56806hx/spins-cv/output
#          (overlays/, per_image/, summary/run_log.txt) + this job's slurm log.
# Pos:     stage-1 (heavy) entrypoint of the two-stage pipeline; run_aggregate
#          is run separately afterwards.
#SBATCH -p multicore
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 02:00:00
#SBATCH -J fibrecv_all
#SBATCH -o /net/scratch/j56806hx/spins-cv/output/summary/sbatch_measure_%j.log

cd /net/scratch/j56806hx/spins-cv/fibrecv
uv run python -m fibrecv.run_measure --all --jobs 8
