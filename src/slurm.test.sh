#!/usr/bin/env bash
#SBATCH --gres=gpu:1
#SBATCH --time=0-00:30:00
#SBATCH --sockets-per-node=1
#SBATCH --partition=slurm_shortgpu
#SBATCH --ntasks=1 --cpus-per-task=30
cd $SLURM_SUBMIT_DIR

# Usage
# ./collide [meshfile] [spherefile] [radius] [outfile]

# Sample Test
MESHFILE=$SLURM_SUBMIT_DIR/../examples/meshfile/sample_mesh.obj
SPHEREFILE=$SLURM_SUBMIT_DIR/../examples/spherefile/sample_spheres.csv
OUTFILE=data.txt
RADIUS=3

./collide $MESHFILE $SPHEREFILE $RADIUS $OUTFILE