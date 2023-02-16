import argparse
import logging
import os
from collections import Counter

# import subprocess
# import sys
import time

# from datetime import datetime, timedelta
from multiprocessing.dummy import Pool
from pathlib import Path
from typing import List

import pandas as pd

from pymatgen.core import Structure
from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor

from vasp_volume_opt import optimize_volume

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__file__)

parser = argparse.ArgumentParser(
    description="Parallel calculation of optimized volumes"
)
parser.add_argument("-i", "--input-structures", type=str, required=True)
#parser.add_argument("-v", "--volume-predictions", type=str, required=True,
#        help="CSV of volume predictions to apply to input structures before running the volume optimization")
parser.add_argument('--vol-pred-site-bias', type=Path,
                    help='Apply a volume prediction to the structures '
                    'before performing the volume optimization. '
                    'Give the path to a file with the average volume per element')
parser.add_argument("-I", "--index", type=int, 
        help="index of structure to run (used for more easily running structures in parallel)")
parser.add_argument("-o", "--out-dir", type=Path,
        default=Path("/projects/rlmolecule/jlaw/volume_relaxation_outputs/")
        )
parser.add_argument("-n", "--poolsize", type=int, required=True)
parser.add_argument("-b", "--batch", type=int, default=0, required=False)


def comptype_to_str(x: List[int]) -> str:
    return "_" + "_".join((str(i) for i in x))


def optimize_row(row: pd.Series) -> None:

    # # Check the time left on the job to avoid submitting a bunch of errors at the end
    # time_left = (
    #     subprocess.check_output('squeue -h -j $SLURM_JOBID -o "%L"', shell=True)
    #     .decode(sys.stdout.encoding)
    #     .strip()
    # )
    # while len(time_left) < 7:
    #     time_left = "00:" + time_left

    # logger.debug(f"Starting {row.id = }, {time_left = }")
    # t = datetime.strptime(time_left,
    #  "%d-%H:%M:%S" if "-" in time_left else "%H:%M:%S")
    # delta = timedelta(hours=t.hour, minutes=t.minute, seconds=t.second)
    # if delta.total_seconds() < 60:
    #     raise RuntimeError("Job Ending")

    try:
        optimize_volume(
            row.unrel_strc_predvol, 
            row.dls_volume, 
            row.id, 
            #row.comptype if '_' in row.comptype else comptype_to_str(row.comptype),
            output_root,
        )

    except Exception as ex:
        logger.warning(f"Exception with {row.id}: {str(ex)}")
        time.sleep(60)


def scale_by_pred_vol(structure: Structure) -> Structure:
    # first predict the volume using the average volume per element (from ICSD)
    site_counts = pd.Series(Counter(
        str(site.specie) for site in structure.sites)).fillna(0)
    curr_site_bias = vol_pred_site_bias[
        vol_pred_site_bias.index.isin(site_counts.index)]
    linear_pred = site_counts @ curr_site_bias
    structure.scale_lattice(linear_pred)

    # then apply Pymatgen's DLS predictor
    pred_volume = dls_vol_predictor.predict(structure)
    structure.scale_lattice(pred_volume)
    return structure


if __name__ == "__main__":
    args = parser.parse_args()
    print(f"Reading {args.input_structures}")
    to_calculate = pd.read_pickle(args.input_structures)
    print(f"\t{len(to_calculate)} structures")

    #to_calculate = to_calculate.reset_index().rename(columns={'index': 'id'})
    output_root = args.out_dir

    # strip previous calculations
    previous_runs = pd.Series(Path(x).name for x in output_root.glob("*/*"))
    #err_runs = pd.Series(Path(x).parent.name for x in output_root.glob("*/*/error"))
    #missing_runs = pd.Series(Path(x).name for x in output_root.glob("*/*")
    #    if not os.path.isfile(f"{x}/CONTCAR")) 

    to_calculate = to_calculate[~to_calculate.id.isin(previous_runs)]
    #to_calculate = to_calculate[to_calculate.id.isin(missing_runs)]
    logger.info(
        f"Restarting {len(to_calculate)}"
        " previous calculations"
    )

    if args.vol_pred_site_bias:
        print("Making linear + DLS volume predictions")
        print(f"Reading {args.vol_pred_site_bias}")
        vol_pred_site_bias = pd.read_csv(args.vol_pred_site_bias,
                                index_col=0, squeeze=True)
        print(f"\t{len(vol_pred_site_bias)} elements")
        dls_vol_predictor = DLSVolumePredictor()
        to_calculate['unrel_strc_predvol'] = to_calculate.structure.apply(scale_by_pred_vol)
        to_calculate['dls_volume'] = to_calculate['unrel_strc_predvol'].apply(lambda x: x.volume)
    if args.index is not None:
        print(args.index)
        to_calculate = to_calculate.iloc[[args.index]]
    else:
        # randomize the order to try to finish more calculations
        to_calculate = to_calculate.sample(frac=1)
    print(to_calculate.head(2))

    # to_calculate = to_calculate[to_calculate.batch == args.batch]
    # logger.info(f"Running batch {args.batch}: {len(to_calculate)}")

    # to_calculate = to_calculate.head(10)

    with Pool(args.poolsize) as p:
        p.map(optimize_row, (row for _, row in to_calculate.iterrows()))
