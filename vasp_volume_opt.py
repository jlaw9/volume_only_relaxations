import csv
import logging
import os
import shutil
import socket
import subprocess
import time
from pathlib import Path
from typing import Tuple
import socket

import yaml
from ase.calculators.vasp import Vasp
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from scipy import optimize

logging.basicConfig(
        format='%(asctime)s:%(levelname)s:%(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

os.environ["VASP_PP_PATH"] = "/projects/rlmolecule/pstjohn/vasp_pp"
vasp_cmd = (
    # "srun --ntasks 18 --cpu-bind rank --exclusive -t 10:00"
    # the -t 10:00 option will limit this subtask to 10 minutes
    #"srun --nodes 1 --ntasks 36 --exclusive -t 10:00"
    "srun --nodes 1 --ntasks 36 --exclusive"
    " /scratch/pgorai/cheng-wei/vasp_std"
)

if os.path.exists(os.environ["LOCAL_SCRATCH"]):
    scratch = os.environ["LOCAL_SCRATCH"] + "/"
else:
    scratch = "/scratch/"


curr_dir = Path(__file__).parent.absolute()
with open(Path(curr_dir, "vasp_opts.yaml"), "r") as f:
    vasp_opts = yaml.safe_load(f)


def calc_energy(volume: float, 
                structure: Structure, 
                directory: str, 
                id_: str,
                file_handle=None,
                ) -> float:

    start = time.time()
    structure.scale_lattice(volume)
    atoms = AseAtomsAdaptor().get_atoms(structure)

    if "F" in id_:
        vasp_opts["encut"] = 540.0

    vasp = Vasp(directory=directory, command=vasp_cmd, **vasp_opts)
    atoms.calc = vasp
    atoms.calc.calculate(atoms)
    energy_zero = atoms.calc.energy_zero

    logger.info(f"SCF ({socket.gethostname()}): {id_ = }, {volume = }, {energy_zero = }")

    f = file_handle
    if file_handle is None:
        f = open(Path(directory, "scf_history.csv"), "a")
    writer = csv.writer(f, delimiter=",")
    writer.writerow([volume, energy_zero, time.time() - start])
    if file_handle is None:
        f.close()

    return energy_zero


def optimize_volume(
    structure: Structure, 
    dls_volume: float, 
    id_: str, 
    #comptype: str, 
    output_root: Path,
) -> Tuple[float, float]:

    # with TemporaryDirectory(prefix=scratch) as tmpdir:
    tmpdir = Path("/scratch/jlaw/20220813_volume_calcs/", id_)
    if tmpdir.exists():
        shutil.rmtree(tmpdir)

    tmpdir.mkdir(parents=True, exist_ok=True)  # exists OK for restarting job

    output_dir = Path(output_root, id_)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:

        # dls_volume = pred_volume(structure)

        assert os.path.exists(tmpdir), f"Issue with tmpdir on {socket.gethostname()}"

        #file_handle = open(Path(tmpdir, "scf_history.csv"), "a")
        file_handle = None

        result = optimize.minimize_scalar(
            calc_energy,
            bounds=(10.0, 2 * dls_volume),
            args=(structure, str(tmpdir), id_, file_handle),
            method="bounded",
            options={"maxiter": 50, "xatol": 1.0},
        )

    except Exception as ex:
        with open(Path(output_dir, "error"), "wt") as f:
            f.write(str(ex))
        return

    #file_handle.close()
    shutil.copy(Path(tmpdir, "scf_history.csv"), output_dir)

    if result.success:
        for file_ in ["INCAR", "CONTCAR", "OUTCAR"]:
            filepath = Path(tmpdir, file_).absolute()
            # leave the CONTCAR file unzipped for now
            if file_ == "CONTCAR":
                shutil.copy(filepath, output_dir)
            else:
                subprocess.run(["gzip", filepath])
                shutil.copy(filepath.with_suffix(".gz"), output_dir)

        with open(Path(output_dir, "opt_vol"), "wt") as f:
            f.write(f"{result.x}\n")
            f.write(f"{result.fun}\n")

    else:
        with open(Path(output_dir, "error"), "wt") as f:
            f.write(f"{result.message}\n")
