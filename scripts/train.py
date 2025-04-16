import numpy as np
from ase.mep.neb import NEB
from ase.optimize import FIRE
from mace.calculators.mace import MACECalculator
from ase.io import read, write
from ase.neighborlist import NeighborList, natural_cutoffs, get_connectivity_matrix
from orca_eval import single_point
import os
import subprocess
import time
import argparse
import logging


logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reactants_fname", help="name of training set file", default="r.xyz")
    parser.add_argument("--products_fname", help="name of training set file", default="p.xyz")
    parser.add_argument("--train_num_rattle", type=int, default=20)
    parser.add_argument("--valid_num_rattle", type=int, default=20)
    parser.add_argument("--test_num_rattle", type=int, default=20)
    parser.add_argument("--rattle_stdev", type=int, default=0.02)
    parser.add_argument("--charge", help="charge of the system", type=int, default=0)
    parser.add_argument("--mult", help="spin multiplicity of the system", type=int, default=1)
    parser.add_argument("--orca_path", help="path top ORCA binary", required=True)
    parser.add_argument("--mace_dir", help="path top MACE installation", required=True)
    parser.add_argument("--num_models", help="number of models in committee (minimum 2)", type=int, default=2)
    parser.add_argument("--delta_threshold", help="model energy disagreement threshold for selection of new data (eV/A)", type=int, default=0.0005)
    parser.add_argument("--num_new_configs", help="max number of configurations added to the training set per iteration", type=int, default=1)
    parser.add_argument("--neb_steps", help="max number of NEB steps", type=int, default=250)
    parser.add_argument("--num_images", help="number of NEB images", type=int, default=20)
    parser.add_argument("--k", help="NEB spring constant", type=float, default=0.2)
    parser.add_argument("--neb_fmax", help="NEB force stopping criterion", type=float, default=0.05)
    parser.add_argument("--qm_fmax", help="Don't accept configurations with force QM force error exceeding this value", type=float, default=10)
    parser.add_argument("--max_n_iterations", help="maximum number of active learning iterations", type=int, default=30)
    parser.add_argument("--extra_neb_steps", action="store_true")
    parser.add_argument("--use_climbing_image", help="Use climbing image NEB after simple NEB", action="store_true")
    parser.add_argument("--use_idpp", help="Use image-dependent pair potential for NEB", action="store_true")
    parser.add_argument("--neb_method", help="which neb method to use: aseneb, improvedtangent", default='aseneb')
    parser.add_argument("--mace_hidden_irreps", help="Number of MACE embedding channels")
    parser.add_argument("--E0s", help="String containing E0s", required=True)
    return parser.parse_args()


def atoms_overlap(atoms):
    cutoffs = natural_cutoffs(atoms, mult=0.6)
    nl = NeighborList(cutoffs, self_interaction=False, skin=0)
    nl.update(atoms)
    return np.any(nl.get_connectivity_matrix(sparse=False))


def initialize_set(fname, r_atoms, p_atoms, indices, stdev=0.02, **kwargs):
    atoms_list = []
    for i in indices:
        atoms = r_atoms.copy()
        atoms.rattle(stdev=stdev, seed=2*i)
        atoms = single_point(config=atoms, **kwargs)
        atoms_list.append(atoms)

        atoms = p_atoms.copy()
        atoms.rattle(stdev=stdev, seed=2*i+1)
        atoms = single_point(config=atoms, **kwargs)
        atoms_list.append(atoms)
    write(fname, atoms_list)


def committee_disagreement(dyn, threshold, neb_steps=250):
    if dyn.nsteps < neb_steps:
        return 0
    images = dyn.atoms.images[1:-1]
    global_threshold = threshold * len(images[0])
    disagreements = np.array([np.sqrt(image.calc.results["energy_var"]) for image in images])
    logger.info(f'Max disagreement: {np.max(disagreements)}, disagreements: {disagreements}')
    logger.info(f'Number of possible new data: {np.sum((disagreements > global_threshold) & ~np.array([atoms_overlap(at) for at in images]))}')
    if np.any((disagreements > global_threshold) & ~np.array([atoms_overlap(at) for at in images])):
        dyn.nsteps = dyn.max_steps


"""
    Selection function that combines mace energy and mace disagreement
    Takes dyn, returns selected images with ref energies and forces
"""


def energy_disagreemenent_selection(dyn, n_select, threshold, orca_path, charge, mult, qm_fmax):
    intermed_images = dyn.atoms.images[1:-1]
    disagreements = np.array([np.sqrt(image.calc.results["energy_var"]) for image in intermed_images])
    images = [image for disagreement, image in zip(disagreements, intermed_images) if disagreement > threshold*len(image)]
    e_score = np.array([np.sqrt(image.calc.results["energy_var"]) for image in images])
    e_score_sorted = np.flip(np.sort(e_score))
    """ List containing indices of images in the order of priority """
    images_order = [np.where(e_score == score)[0][0] for score in e_score_sorted]
    selected_images = []
    for index in images_order:
        if len(selected_images) == n_select:
            break
        image = images[index]
        
        if not atoms_overlap(image):
            image = single_point(orca_path=orca_path, config=image, charge=charge, mult=mult)
        else:
            continue

        if np.all(np.linalg.norm(image.arrays["mp2_forces"], axis=1) <= qm_fmax):
            logger.info('Selected image')
            selected_images.append(image)
        else:
            logger.info('QM forces too high - going to next image')
            continue

    return selected_images


def wait_for_output(fnames, sleep_int=10):
    while True:
        if all([os.path.exists(fname) for fname in fnames]):
            break
        time.sleep(sleep_int)


def initialize_neb(reac, prod, nimages=10, k=0.1, neb_method="aseneb", use_idpp=False):
    images = ([reac] 
            + [reac.copy() for i in range(nimages-2)] 
            + [prod])
    neb = NEB(images=images, k=k, remove_rotation_and_translation=True, method=neb_method)
    if use_idpp:
        neb.interpolate('idpp')
    else:
        neb.interpolate()
    return neb


def write_neb(dyn):
    for im in dyn.atoms.images:
        if 'energies' in im.calc.results:
            del im.calc.results['energies']
        if 'forces_comm' in im.calc.results:
            del im.calc.results['forces_comm']
    write('neb.xyz', dyn.atoms.images)


def main():
    logging.basicConfig(filename='al.log', level=logging.INFO)
    logger.info('Started')
    
    args = parse_args()
    
    reac = read(args.reactants_fname)
    prod = read(args.products_fname)
    mace_dir = args.mace_dir
    orca_path = args.orca_path
    num_models = args.num_models
    n_select = args.num_new_configs
    neb_steps = args.neb_steps
    num_images = args.num_images
    k = args.k
    neb_fmax = args.neb_fmax
    qm_fmax = args.qm_fmax
    max_n_iterations = args.max_n_iterations
    charge = args.charge
    mult = args.mult
    model = "MACE_*_stagetwo.model"

    logger.info(f'Using {num_images} images in NEB')
    logger.info(f'Using spring constant of {k}')
    logger.info(f'Using stopping force threshold of {neb_fmax} eV / A and maximum of {neb_steps} steps')
    logger.info(f'Using maximum true force of {qm_fmax} for accepting new configs')
    logger.info(f'Taking maximum of {n_select} new configurations per iteration')

    # Check for train, valid, test sets. If not present, build them
    train_file_exists = os.path.exists('train.xyz')
    valid_file_exists = os.path.exists('valid.xyz')
    test_file_exists = os.path.exists('test.xyz')

    if train_file_exists:
        logger.info('Found training file train.xyz')
    else:
        logger.info(f'Training file not found. Creating using {args.train_num_rattle} rattled reactant and product configs.')
        indices = np.arange(args.train_num_rattle)
        initialize_set('train.xyz', reac, prod, indices, 
            stdev=args.rattle_stdev, orca_path=orca_path, charge=charge, mult=mult
        )
    
    if valid_file_exists:
        logger.info('Found validation file valid.xyz')
    else:
        logger.info(f'Validation file not found. Creating using {args.valid_num_rattle} rattled reactant and product configs.')
        indices = np.arange(args.train_num_rattle, args.train_num_rattle+args.valid_num_rattle)
        initialize_set('valid.xyz', reac, prod, indices, 
            stdev=args.rattle_stdev, orca_path=orca_path, charge=charge, mult=mult
        )
    
    if test_file_exists:
        logger.info('Found test file test.xyz')
    else:
        logger.info(f'Test file not found. Creating using {args.test_num_rattle} rattled reactant and product configs.')
        indices = np.arange(args.train_num_rattle+args.valid_num_rattle, args.train_num_rattle+args.valid_num_rattle+args.test_num_rattle)
        initialize_set('test.xyz', reac, prod, indices, 
            stdev=args.rattle_stdev, orca_path=orca_path, charge=charge, mult=mult
        )

    for iter in range(max_n_iterations):
        if os.path.exists(f"iter{iter}"):
            logger.info(f"Iteration {iter} already done. Skipping to next iteration.")
            continue
        commands = []
        for s in range(1, num_models+1):
            if os.path.exists(f"MACE_{s}_stagetwo.model"):
                logger.info(f"{s}th model already trained for current iteration. Skipping to next seed.")
                continue
            seed = iter * num_models + s
            commands.append([
                "python", f"{mace_dir}/cli/run_train.py",
                f"--name=MACE_{s}", "--train_file=train.xyz", 
                "--valid_file=valid.xyz", "--test_file=test.xyz",
                "--model=MACE", "--loss=weighted", f"--E0s={args.E0s}",
                f"--hidden_irreps={args.mace_hidden_irreps}",
                "--batch_size=5", "--swa", "--ema",
                "--max_num_epochs=400", "--start_swa=200",
                "--ema_decay=0.99", "--amsgrad",
                "--restart_latest", "--device=cuda", 
                "--save_cpu", f"--seed={seed}",
                "--forces_key=REF_forces", "--energy_key=REF_energy",
                "--default_dtype=float64"])
        if len(commands) == 0:
            logger.info("All models available, continuing to NEB")
        else:
            processes = [subprocess.Popen(command) for command in commands]
            for process in processes:
                process.wait()
        wait_for_output([f"MACE_{s}_stagetwo.model" for s in range(1, num_models+1)])
        logger.info('Doing NEB with MACE models to get new training data')
 
        if not os.path.exists("neb.xyz"):
            neb = initialize_neb(reac=reac, prod=prod, nimages=num_images, k=k, neb_method=args.neb_method, use_idpp=args.use_idpp)
        else:
            neb = NEB(images=read("neb.xyz", ":"), k=k, method='improvedtangent', remove_rotation_and_translation=True)
    
        # Set calculators:
        for image in neb.images:
            image.calc = MACECalculator(model_paths=model, default_dtype="float64", device="cuda")
        # Run neb for neb_steps initial steps:
        dyn = FIRE(neb)
        dyn.run(fmax=neb_fmax, steps=args.neb_steps)
        dyn.attach(write_neb, interval=1, dyn=dyn)

        if args.extra_neb_steps:
            dyn.attach(committee_disagreement, interval=1, dyn=dyn, threshold=0.0005, neb_steps=args.neb_steps)
            dyn.run(fmax=neb_fmax, steps=neb_steps)
        
        if args.use_climb:
            dyn.atoms.climb = True
            dyn.run(fmax=neb_fmax, steps=dyn.nsteps + args.neb_steps)

        """ Select new config """
        selected_images = energy_disagreemenent_selection(dyn=dyn, n_select=n_select, threshold=0.0005, orca_path=orca_path, charge=charge, mult=mult, qm_fmax=qm_fmax)
        logger.info(f"number of selected configs {len(selected_images)}")

        """ Create new directory for model, neb.xyz and copy of train and valid """
        os.system(f"mkdir iter{iter}")
        os.system("rm MACE*compiled.model")
        os.system(f"mv neb.xyz MACE* logs iter{iter}")
        os.system(f"cp train.xyz iter{iter}")
        os.system("rm -r results checkpoints")
        

        if len(selected_images) > 0:
            for im in selected_images:
                if 'energies' in im.calc.results.keys():
                    del im.calc.results['energies']
            write("train.xyz", selected_images, append=True)
    
        else:
            logger.info("No new configuration selected. AL has converged.")
            break


if __name__ == "__main__":
    main()
