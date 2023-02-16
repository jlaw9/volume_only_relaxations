#!/bin/bash
#SBATCH --account=rlmolecule
#SBATCH --time=1:00:00
#SBATCH --job-name=test_ase_async
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --partition=debug

source ~/.bashrc
module load gcc mkl
conda activate /home/jlaw/.conda-envs/crystals_nfp0_3

#python run_volume_optimization.py -n 2 -b 1
#python -u run_volume_optimization.py \
#       --input-structures=/projects/rlmolecule/jlaw/crystal_outputs/2022-07-05/volume_relaxation_outputs/best_decorations_volpred.p \
#       --out-dir=/projects/rlmolecule/jlaw/crystal_outputs/2022-07-05/volume_relaxation_outputs \
#       --poolsize=1
python -u run_volume_optimization.py \
    --input-structures=/projects/rlmolecule/jlaw/inputs/structures/battery/all_unrelaxed_structures.p \
    --vol-pred-site-bias=/projects/rlmolecule/pstjohn/crystal_inputs/site_volumes_from_icsd.csv \
    --poolsize=1
