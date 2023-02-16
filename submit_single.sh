#!/bin/bash
#SBATCH --account=rlmolecule
##SBATCH --time=2-00
#SBATCH --time=4:00:00
#SBATCH --job-name=ase_multiprocess
#SBATCH --nodes=25
#SBATCH --ntasks-per-node=36
#SBATCH --mail-user=jlaw@nrel.gov
#SBATCH --mail-type=BEGIN,END

source ~/.bashrc
module load gcc mkl
conda activate /home/jlaw/.conda-envs/crystals_nfp0_3

echo "Job started: `date`"
python -u run_volume_optimization.py \
       --input-structures=/projects/rlmolecule/jlaw/crystal_outputs/2022-07-05/volume_relaxation_outputs/best_decorations_volpred.p \
       --out-dir=/projects/rlmolecule/jlaw/crystal_outputs/2022-07-05/volume_relaxation_outputs \
       --poolsize=25

echo "Job finished: `date`"
