import os
from use_clean import activate_alphagem_env
import subprocess
import shutil
import sys

def activate_deepectransformer_env():
    try:
        subprocess.run("source activate deepectransformer", shell=True, check=True)
        print("Switched to deepectransformer environment.")
    except subprocess.CalledProcessError as e:
        print(f"Error switching to CLEAN environment: {e}")


def use_deepectransformer(name):
    path = os.getcwd()
    original_cwd = path
    tools_dir = os.path.join(path, "tools", "DeepProZyme")
    script_path = os.path.join(path, "tools", "DeepProZyme", "run_deepectransformer.py")
    input_fasta = os.path.join(path, "working", name, f"{name}_homoleft.fasta")
    output_dir = os.path.join(path, "working", name, f"{name}_deepec_result")

    if not os.path.isdir(tools_dir):
        print(f"DeepProZyme directory not found: {tools_dir}")
        return
    if not os.path.isfile(input_fasta):
        print(f"Input FASTA not found: {input_fasta}")
        return

    os.makedirs(output_dir, exist_ok=True)

    def _exists(cmd: str) -> bool:
        return shutil.which(cmd) is not None

    def _candidate_cmds():
        cmds = []
        conda_exe = os.environ.get("CONDA_EXE")
        if conda_exe and os.path.isfile(conda_exe):
            cmds.append([conda_exe, "run", "-n", "deepectransformer", "python"])
        if _exists("conda"):
            cmds.append(["conda", "run", "-n", "deepectransformer", "python"])
        if _exists("mamba"):
            cmds.append(["mamba", "run", "-n", "deepectransformer", "python"])
        if _exists("micromamba"):
            cmds.append(["micromamba", "run", "-n", "deepectransformer", "python"])
        cmds.append([sys.executable])
        return cmds

    try:
        os.chdir(tools_dir)
        #activate_deepectransformer_env()
        base_args = [
            script_path,
            "-i", input_fasta,
            "-o", output_dir,
            "-g", "cpu",
            "-b", "128",
            "-cpu", "2",
        ]
        last_err = None
        for runner in _candidate_cmds():
            cmd = runner + base_args
            try:
                print("Running:", " ".join(cmd))
                subprocess.run(cmd, check=True)
                print("deepectransformer completed.")
                break
            except subprocess.CalledProcessError as e:
                last_err = e
                print(f"Runner failed: {' '.join(runner)}; return code {e.returncode}")
        else:
            if last_err:
                raise last_err
    finally:
        os.chdir(original_cwd)
        #activate_alphagem_env()