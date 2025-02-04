from ase.io import read,write
import os
import subprocess
import numpy as np
import time
import glob
import re
import argparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_fname", help="name of training set file", default="train.xyz")
    parser.add_argument("--test_fname", help="name of training set file", default="test.xyz")
    parser.add_argument("--train_eval_fname", help="file to evaluate training set MACE energies", default="train_eval.xyz")
    parser.add_argument("--test_eval_fname", help="file to evaluate test set MACE energies", default="test_eval.xyz")
    parser.add_argument("--mace_dir", help="path for MACE installation", required=True)
    parser.add_argument("--num_new_configs", help="number of configurations to transfer to the training set per iteration", type=int, default=45)
    parser.add_argument("--num_iterations", help="number of iterations", type=int, default=15)
    return parser.parse_args()

def wait_for_output(fnames, sleep_int=10):
    output_not_here_yet = True
    while output_not_here_yet:
        all_outputs_exist = True
        for fname in fnames:
            all_outputs_exist = all_outputs_exist and os.path.exists(fname)
        if all_outputs_exist:
            output_not_here_yet = False
        time.sleep(sleep_int)
        

def update_set(test_fname, num_new_configs):
    test = read(test_fname, ":")
    energy_error = np.array([np.abs(config.info["mp2_energy"] - config.info["MACE_energy"]) for config in test])
    selected_configs = [a for a, b in zip(test, list(energy_error >= np.sort(energy_error)[-num_new_configs])) if b]
    untouched_configs = [a for a, b in zip(test, list(energy_error < np.sort(energy_error)[-num_new_configs])) if b]
    return selected_configs, untouched_configs

def main():
    args = parse_args()
    ##define input variables:
    train_fname = args.train_fname #will contain a few reactant and ts structures
    test_fname = args.test_fname #will contain reactant and ts structures
    test_eval_fname = args.test_eval_fname
    mace_dir = args.mace_dir
    num_new_configs = args.num_new_configs
    # define number of iterations
    num_iterations = args.num_iterations

    for i in range(num_iterations):
        if os.path.exists(f'iter{i}'):
            continue
        folder_path = "iter{}".format(i)
        if not os.path.exists(folder_path):
            subprocess.run(["mkdir", folder_path])

        if i == 0:
            max_num_epochs = 1000
        elif i > 0:
            print(glob.glob(f"checkpoints/MACE_run-{i}_epoch-*_swa.pt"))
            max_num_epochs = 200 + int(re.split("_|-", glob.glob(f"checkpoints/MACE_run-0_epoch-*_swa.pt")[0])[4])
        
        if not all([os.path.exists(folder_path+"/"+fname) for fname in ["MACE.model", "MACE_swa.model"]]):
            if not all([os.path.exists(fname) for fname in ["MACE.model", "MACE_swa.model"]]):
                # submit slurm job to train model
                subprocess.run([
                    "python3", f"{mace_dir}/mace/cli/run_train.py",
                    f"--name=MACE", "--train_file=train.xyz", 
                    "--valid_fraction=0.1", "--test_file=test.xyz",
                    "--model=MACE", "--loss=weighted",
                    "--batch_size=10", "--swa", "--ema",
                    f"--max_num_epochs={max_num_epochs}", "--start_swa=500",
                    "--ema_decay=0.99", "--amsgrad",
                    "--restart_latest", "--device=cuda", 
                    "--save_cpu", f"--seed=0",
                    "--forces_key=mp2_forces", "--energy_key=mp2_energy",
                    "--default_dtype=float64"
                ])
                # wait for the train outputs to appear
                wait_for_output(fnames=["MACE.model", "MACE_swa.model"])
        # evaluate test (preferably save the xyz evaluated at the end of the training)
        if not all([os.path.exists(folder_path+"/"+fname) for fname in ["train_eval.xyz", "test_eval.xyz"]]):
            if not all([os.path.exists(fname) for fname in ["train_eval.xyz", "test_eval.xyz"]]):
                subprocess.run([
                    "python3", f"{mace_dir}/mace/cli/eval_configs.py",
                    "--configs=train.xyz", "--output=train_eval.xyz", 
                    "--model=MACE_swa.model", "--device=cuda", 
                    "--default_dtype=float64"
                ])
                subprocess.run([
                    "python3", f"{mace_dir}/mace/cli/eval_configs.py",
                    "--configs=test.xyz", "--output=test_eval.xyz", 
                    "--model=MACE_swa.model", "--device=cuda", 
                    "--default_dtype=float64"
                ])

                wait_for_output(fnames=["train_eval.xyz", "test_eval.xyz"])


        if i != num_iterations-1:
            # select new configs, update training set, delete old model, update slurm file to set the correct number of epochs
            train = read(train_fname, ":")
            selected_configs, untouched_configs = update_set(test_eval_fname, num_new_configs)
            new_train = train + selected_configs
            write(train_fname, new_train)
            write(test_fname, untouched_configs)
            # number of epochs from the last checkpoint
        
        subprocess.run(["mv", "logs", "results", "MACE.model", "MACE_swa.model", "train_eval.xyz", "test_eval.xyz", folder_path])
 

if __name__ == "__main__":
    main()
