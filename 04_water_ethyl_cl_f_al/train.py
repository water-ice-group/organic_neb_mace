import numpy as np
from ase.mep.neb import NEB, NEBOptimizer
from ase.optimize import FIRE
import glob
from mace.calculators.mace import MACECalculator
from ase.io import read, write
from orca_eval import single_point
import ase
import os
import subprocess
import time
import re
import argparse

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reactants_fname", help="name of training set file", default="r.xyz")
    parser.add_argument("--products_fname", help="name of training set file", default="p.xyz")
    parser.add_argument("--mace_dir", help="path top MACE installation", required=True)
    parser.add_argument("--orca_path", help="path top ORCA binary", required=True)
    parser.add_argument("--num_models", help="number of models in committee (minimum 2)", type=int, default=2)
    parser.add_argument("--num_new_configs", help="max number of configurations added to the training set per iteration", type=int, default=1)
    parser.add_argument("--neb_steps", help="number of NEB steps", type=int, default=200)
    parser.add_argument("--n_images", help="number of NEB images", type=int, default=20)
    parser.add_argument("--k", help="NEB spring constant", type=float, default=0.02)
    parser.add_argument("--neb_fmax", help="NEB force stopping criterion", type=float, default=0.05)
    parser.add_argument("--qm_fmax", help="Configurations exceeding QC force cutoff (in eV / A) are not included in the training set", type=float, default=10)
    parser.add_argument("--max_n_iterations", help="maximum number of active learning iterations", type=int, default=20)
    parser.add_argument("--charge", help="charge of the system", type=int, default=0)
    parser.add_argument("--mult", help="spin multiplicity of the system", type=int, default=1)
    return parser.parse_args()


def committee_error(dyn):
    images = dyn.atoms.images[1:-1]
    errors = []
    errors = [np.sqrt(image.calc.results["energy_var"]) for image in images]
    dyn.atoms.errors = np.array(errors)

"""
    Selection function that combines mace energy and mace disagreement
    Takes dyn, returns selected images with ref energies and forces
"""

def energy_disagreemenent_selection(dyn, n_select, threshold, orca_path, charge, mult):
    images = [image for disagreement, image in zip(dyn.atoms.errors, dyn.atoms.images[1:-1]) if disagreement > threshold*image.get_global_number_of_atoms()]
    e_score = np.array([np.sqrt(image.calc.results["energy_var"]) for image in images])
    e_score_sorted = np.flip(np.sort(e_score))
    """ List containing indices of images in the order of priority """
    images_order = [np.where(e_score == score)[0][0] for score in e_score_sorted]
    selected_images = []
    for index in images_order:
        if len(selected_images) == n_select:
            break
        image = images[index]
        image = single_point(orca_path=orca_path, config=image, charge=charge, mult=mult)
        if np.any(np.linalg.norm(image.arrays["mp2_forces"], axis=1) <= 10):
            selected_images.append(image)
        else:
            continue
    return selected_images

def wait_for_output(fnames, sleep_int=10):
    while True:
        if all([os.path.exists(fname) for fname in fnames]):
            break
        time.sleep(sleep_int)

def initialize_neb(reac, prod, nimages=10, k=0.1):
    images = ([reac] 
            + [reac.copy() for i in range(nimages-2)] 
            + [prod])
    neb = NEB(images=images, k=k, climb=False, remove_rotation_and_translation=True, parallel=True)
    neb.interpolate()
    return neb


def main():
    args = parse_args()
    
    reac = read(args.reactants_fname)
    prod = read(args.products_fname)
    mace_dir = args.mace_dir
    orca_path = args.orca_path
    n_comm = args.num_models
    n_select = args.num_new_configs
    neb_steps = args.neb_steps
    n_images = args.n_images
    k = args.k
    neb_fmax = args.neb_fmax
    qm_fmax = args.qm_fmax
    max_n_iterations = args.max_n_iterations
    charge = args.charge
    mult = args.mult
    model = "MACE_*_swa.model"

    

    for iter in range(max_n_iterations):
        if os.path.exists(f"iter{iter}"):
            print(f"Iteration {iter} aklready done. Skipping to next iteration.")
            continue
        commands = []
        for s in range(1, n_comm+1):
            if os.path.exists(f"MACE_{s}_swa.model"):
                print(f"{s}th model already trained for current iteration. Skipping to next seed.")
                continue
            seed = iter * n_comm + s
            commands.append([
                "python", f"{mace_dir}/mace/cli/run_train.py",
                f"--name=MACE_{s}", "--train_file=train.xyz", 
                "--valid_file=valid.xyz", "--test_file=test.xyz",
                "--model=MACE", "--loss=weighted",
                "--batch_size=5", "--swa", "--ema",
                "--max_num_epochs=1000", "--start_swa=500",
                "--ema_decay=0.99", "--amsgrad",
                "--restart_latest", "--device=cuda", 
                "--save_cpu", f"--seed={seed}",
                "--forces_key=mp2_forces", "--energy_key=mp2_energy",
                "--default_dtype=float64"])
        if len(commands) == 0:
            print("All models available, continuing to NEB")
        else:
            processes = [subprocess.Popen(command) for command in commands]
            for process in processes:
                process.wait()
        wait_for_output([f"MACE_{s}_swa.model" for s in range(1, n_comm+1)])
        print("Do neb")
        if not os.path.exists("neb.xyz"):
            neb = initialize_neb(reac=reac, prod=prod, nimages=n_images, k=k)
        else:
            neb = NEB(images=read("neb.xyz", ":"), k=k)
    
        # Set calculators:
        for image in neb.images[1:-1]:
            image.calc = MACECalculator(model_paths=model, default_dtype="float64", device="cuda")
        # Optimize:
        dyn = FIRE(neb)
        dyn.attach(committee_error, interval=1, dyn=dyn)
        """ Run NEB """
        dyn.run(fmax=neb_fmax, steps=neb_steps)
        write("neb.xyz", dyn.atoms.images)
        """ Select new config """
        selected_images = energy_disagreemenent_selection(dyn=dyn, n_select=n_select, threshold=0.0005)

        """ Create new directory for model, neb.xyz and copy of train and valid """
        os.system(f"mkdir iter{iter}")
        os.system(f"mv neb.xyz MACE* logs iter{iter}")
        os.system(f"cp train.xyz iter{iter}")
        os.system("rm -r results checkpoints")
        print(f"number of selected configs {len(selected_images)}")

        if len(selected_images) > 0:
            write("train.xyz", selected_images, append=True)
    
        else:
            print("No new configuration selected. AL has converged.")
            break

if __name__ == "__main__":
    main()