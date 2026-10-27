"""Select one licensed LP/MIP solver before constructing any COBRA models."""

import argparse
import importlib
import json
import os
import subprocess
import sys

import cobra


SOLVER_CHOICES = ("auto", "gurobi", "copt")
_selected_solver = None


def add_solver_arguments(parser):
    parser.add_argument(
        "--solver", choices=SOLVER_CHOICES, default="auto",
        help="LP/MIP solver: auto tries Gurobi then COPT (default: auto)",
    )


def positive_number(value):
    number = float(value)
    if not 0 < number < float("inf"):
        raise argparse.ArgumentTypeError("must be a finite positive number")
    return number


def nonnegative_number(value):
    number = float(value)
    if not 0 <= number < float("inf"):
        raise argparse.ArgumentTypeError("must be a finite nonnegative number")
    return number


def _interface(name):
    if name == "copt":
        from optlang_copt import register_with_cobra
        # Direct probe execution also needs the repository on sys.path.
        from pathlib import Path

        root = str(Path(__file__).resolve().parents[1])
        if root not in sys.path:
            sys.path.insert(0, root)
        from src.copt_compat import apply_copt_compatibility

        interface = register_with_cobra()
        apply_copt_compatibility(interface)
        return interface
    return importlib.import_module("optlang.gurobi_interface")


def _probe(name):
    """Check LP, binary optimization and production-size license in isolation.

    Native license failures can terminate Python, so probes use a subprocess.
    The >2000-variable check rejects size-limited trial licenses early.
    """
    interface = _interface(name)
    for integer in (False, True):
        model = interface.Model()
        model.configuration.verbosity = 0
        model.configuration.timeout = 20
        variables = [interface.Variable(f"x{i}", lb=0, ub=1) for i in range(2101)]
        indicator = interface.Variable("indicator", type="binary" if integer else "continuous", lb=0, ub=1)
        model.add(variables + [indicator])
        model.add(interface.Constraint(sum(variables) - 2101 * indicator, ub=0))
        model.objective = interface.Objective(sum(variables) - indicator, direction="max")
        status = model.optimize()
        if status != "optimal" or abs(model.objective.value - 2100) > 1e-5:
            raise RuntimeError(f"{'MIP' if integer else 'LP'} probe failed: {status}")


def select_solver(requested="auto"):
    """Set COBRA's process-wide default and the inherited subprocess setting."""
    global _selected_solver
    if requested not in SOLVER_CHOICES:
        raise ValueError(f"Unknown solver: {requested}")
    inherited = os.environ.get("ALPHAGEM_SOLVER")
    if requested == "auto" and inherited:
        if inherited not in SOLVER_CHOICES:
            raise ValueError("ALPHAGEM_SOLVER must be auto, gurobi or copt")
        requested = inherited
    candidates = ("gurobi", "copt") if requested == "auto" else (requested,)
    failures = []
    for name in candidates:
        try:
            result = subprocess.run(
                [sys.executable, os.path.abspath(__file__), "--probe", name],
                capture_output=True, text=True, timeout=60,
            )
            if result.returncode:
                detail = (result.stderr or result.stdout).strip().splitlines()
                raise RuntimeError(detail[-1] if detail else f"probe exited {result.returncode}")
            interface = _interface(name)
            cobra.Configuration().solver = interface
        except (Exception, SystemExit) as error:
            failures.append(f"{name}: {error}")
            print(f"AlphaGEM solver unavailable: {failures[-1]}", flush=True)
            continue
        _selected_solver = name
        os.environ["ALPHAGEM_SOLVER"] = name
        print(f"AlphaGEM solver selected: {name}", flush=True)
        return name
    raise RuntimeError(
        "No usable requested LP/MIP solver. Install/license Gurobi or "
        "optlang-copt with COPT. " + "; ".join(failures)
    )


def ensure_solver():
    """Keep standalone pipeline modules and new/copied models on one backend."""
    if _selected_solver is None:
        return select_solver()
    cobra.Configuration().solver = _interface(_selected_solver)
    return _selected_solver


def write_solver_manifest(path):
    from importlib.metadata import version
    from pathlib import Path

    name = ensure_solver()
    packages = ("cobra", "optlang", "optlang-copt", "coptpy") if name == "copt" else ("cobra", "optlang", "gurobipy")
    Path(path).write_text(json.dumps({
        "solver": name, "versions": {package: version(package) for package in packages},
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] != "--probe":
        raise SystemExit("usage: solver_config.py --probe {gurobi,copt}")
    _probe(sys.argv[2])
