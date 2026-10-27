"""Regression tests for process-wide selection and both gapfilling backends."""

import argparse
from importlib.metadata import version
import math
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
from unittest.mock import Mock, PropertyMock, patch

import cobra
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from src import solver_config
from src.gapfilling import (
    configure_gapfiller,
    fill_with_validated_incumbent,
    fill_with_validation_retry,
    has_growth,
    require_growth,
)


@pytest.fixture(autouse=True)
def restore_configuration(monkeypatch):
    original = cobra.Configuration().solver
    monkeypatch.delenv("ALPHAGEM_SOLVER", raising=False)
    monkeypatch.setattr(solver_config, "_selected_solver", None)
    yield
    cobra.Configuration().solver = original


def test_auto_falls_back_after_native_failure():
    results = [SimpleNamespace(returncode=-11, stderr="", stdout=""),
               SimpleNamespace(returncode=0, stderr="", stdout="")]
    with patch.object(solver_config.subprocess, "run", side_effect=results) as run:
        assert solver_config.select_solver() == "copt"
        assert [call.args[0][-1] for call in run.call_args_list] == ["gurobi", "copt"]
    assert os.environ["ALPHAGEM_SOLVER"] == "copt"
    assert cobra.Model().solver.interface.__name__ == "optlang_copt.interface"
    assert cobra.Model().copy().solver.interface.__name__ == "optlang_copt.interface"


def test_explicit_failure_does_not_switch_solver():
    with patch.object(solver_config.subprocess, "run", return_value=SimpleNamespace(
            returncode=1, stderr="license invalid", stdout="")) as run:
        with pytest.raises(RuntimeError, match="license invalid"):
            solver_config.select_solver("gurobi")
    assert run.call_count == 1


def test_auto_respects_environment(monkeypatch):
    monkeypatch.setenv("ALPHAGEM_SOLVER", "copt")
    with patch.object(solver_config.subprocess, "run", return_value=SimpleNamespace(
            returncode=0, stderr="", stdout="")) as run:
        assert solver_config.select_solver() == "copt"
    assert run.call_args.args[0][-1] == "copt"


def test_cli_overrides_environment(monkeypatch):
    monkeypatch.setenv("ALPHAGEM_SOLVER", "gurobi")
    with patch.object(solver_config.subprocess, "run", return_value=SimpleNamespace(
            returncode=0, stderr="", stdout="")) as run:
        assert solver_config.select_solver("copt") == "copt"
    assert run.call_count == 1


def test_copt_021_uses_native_fixes_without_application_patch():
    if version("optlang-copt") != "0.2.1":
        pytest.skip("Native compatibility regressions require optlang-copt 0.2.1")
    interface = solver_config._interface("copt")
    assert not getattr(interface.Model, "_alphagem_rollback_compatibility", False)
    assert hasattr(interface.Model, "_has_primal_solution")


@pytest.mark.parametrize("value", ["0", 0, "0.01"])
def test_nonnegative_number_accepts_zero_and_positive_values(value):
    assert solver_config.nonnegative_number(value) == float(value)


@pytest.mark.parametrize("value", ["-0.01", "nan", "inf", "-inf"])
def test_nonnegative_number_rejects_negative_and_nonfinite_values(value):
    with pytest.raises(argparse.ArgumentTypeError, match="finite nonnegative"):
        solver_config.nonnegative_number(value)


@pytest.mark.parametrize("variable_type", ["continuous", "binary"])
def test_copt_linear_variable_removal_preserves_edited_rows_and_limits(variable_type):
    interface = solver_config._interface("copt")
    model = interface.Model("edited_rows")
    retained = interface.Variable("retained", lb=0, ub=10)
    removed = interface.Variable("removed", lb=0, ub=1, type=variable_type)
    model.add([retained, removed])
    row = interface.Constraint(retained, ub=4, name="capacity")
    model.add(row)
    model.objective = interface.Objective(retained, direction="max")
    # COBRA edits native coefficients after adding reactions to existing rows.
    # The removed variable was absent from the original symbolic constraint.
    row.set_linear_coefficients({retained: 2, removed: 3})
    model.objective.set_linear_coefficients({retained: 2, removed: 5})
    native_problem = model.problem
    native_problem.setParam("RelGap", 0.002)
    native_problem.setParam("MemLimit", 2048)
    assert model.optimize() == "optimal"
    assert model.objective.value == pytest.approx(6)

    model.remove(removed)
    model.update()

    # A fresh native problem would discard these resource settings and make
    # large batches of gene/reaction deletions prohibitively expensive.
    assert model.problem is native_problem
    assert model.problem.getParam("RelGap") == pytest.approx(0.002)
    assert model.problem.getParam("MemLimit") == 2048
    assert list(model.variables) == [retained]
    assert len(model.problem.getVars()) == 1
    assert row.expression.free_symbols == {retained}
    assert model.objective.expression.free_symbols == {retained}
    assert row.get_linear_coefficients([retained])[retained] == 2
    assert model.objective.get_linear_coefficients([retained])[retained] == 2
    assert model.optimize() == "optimal"
    assert retained.primal == pytest.approx(2)
    assert model.objective.value == pytest.approx(4)


def test_copt_quadratic_moma_rolls_back_after_restoring_linear_objective():
    interface = solver_config._interface("copt")
    cobra.Configuration().solver = interface
    model = cobra.Model("quadratic_rollback")
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    uptake = cobra.Reaction("UPTAKE", upper_bound=10)
    uptake.add_metabolites({a: 1})
    conversion = cobra.Reaction("CONVERT")
    conversion.add_metabolites({a: -1, b: 1})
    biomass = cobra.Reaction("BIOMASS")
    biomass.add_metabolites({b: -1})
    model.add_reactions([uptake, conversion, biomass])
    model.objective = biomass
    reference = model.optimize()
    assert reference.objective_value == pytest.approx(10)
    original_variables = set(model.solver.variables.keys())
    original_constraints = set(model.solver.constraints.keys())
    native_problem = model.solver.problem

    with model:
        conversion.upper_bound = 5
        cobra.flux_analysis.add_moma(model, solution=reference, linear=False)
        adjusted = model.optimize()
        assert adjusted.status == "optimal"
        assert adjusted.fluxes["BIOMASS"] == pytest.approx(5)
        assert adjusted.objective_value == pytest.approx(75)
        assert model.solver.problem.HasQObj
        # COBRA's context undo also restores an LP objective before deleting
        # MOMA auxiliaries. Former QP variables still require a safe rebuild.
        model.objective = interface.Objective(biomass.flux_expression, direction="max")
        assert model.solver.objective.is_Linear
        assert not model.solver.problem.HasQObj
        assert model.slim_optimize() == pytest.approx(5)

    assert model.solver.problem is not native_problem
    assert not model.solver.problem.HasQObj
    assert set(model.solver.variables.keys()) == original_variables
    assert set(model.solver.constraints.keys()) == original_constraints
    assert model.slim_optimize() == pytest.approx(10)
    assert model.solver.objective.direction == "max"
    assert model.metabolites.a_c.constraint.get_linear_coefficients(
        [uptake.forward_variable, conversion.forward_variable]
    ) == {uptake.forward_variable: 1, conversion.forward_variable: -1}


def test_copt_restoring_linear_objective_clears_native_quadratic_terms():
    interface = solver_config._interface("copt")
    model = interface.Model("replace_quadratic_objective")
    variable = interface.Variable("x", lb=-1, ub=1)
    model.add(variable)
    model.objective = interface.Objective(variable**2, direction="min")
    assert model.optimize() == "optimal"
    assert model.objective.value == pytest.approx(0, abs=1e-7)

    model.objective = interface.Objective(variable, direction="max")

    assert model.objective.is_Linear
    assert not model.problem.HasQObj
    assert model.optimize() == "optimal"
    assert variable.primal == pytest.approx(1)
    assert model.objective.value == pytest.approx(1)


def test_copt_model_copy_preserves_solver_configuration():
    interface = solver_config._interface("copt")
    cobra.Configuration().solver = interface
    model = cobra.Model("copy_configuration")
    model.tolerance = 1e-9
    model.solver.configuration.timeout = 12
    model.solver.configuration.presolve = False
    model.solver.configuration.lp_method = "dual"

    copied = model.copy()

    assert copied.tolerance == model.tolerance == 1e-9
    assert copied.solver.configuration.tolerances.feasibility == 1e-9
    assert copied.solver.configuration.tolerances.integrality == 1e-9
    assert copied.solver.configuration.tolerances.optimality == 1e-9
    assert copied.solver.configuration.timeout == 12
    assert copied.solver.configuration.presolve is False
    assert copied.solver.configuration.lp_method == "dual"


def test_copt_infeasible_model_returns_missing_fluxes_and_preserves_status():
    interface = solver_config._interface("copt")
    cobra.Configuration().solver = interface
    model = cobra.Model("infeasible_draft")
    metabolite = cobra.Metabolite("a_c", compartment="c")
    production = cobra.Reaction("PRODUCTION", lower_bound=1, upper_bound=1)
    production.add_metabolites({metabolite: 1})
    model.add_reactions([production])
    model.objective = production

    with pytest.warns(UserWarning, match="infeasible"):
        solution = model.optimize()

    assert solution.status == "infeasible"
    assert model.solver.status == "infeasible"
    assert math.isnan(solution.objective_value)
    assert solution.fluxes.index.tolist() == ["PRODUCTION"]
    assert solution.fluxes.isna().all()
    assert not has_growth(model, SimpleNamespace(objective_value=1e30))
    with pytest.raises(cobra.exceptions.OptimizationError, match="infeasible"):
        model.optimize(raise_error=True)
    with pytest.warns(UserWarning, match="infeasible"):
        with pytest.raises(RuntimeError, match="no finite positive growth"):
            require_growth(model, "Test draft")


@pytest.fixture
def solved_copt_gapfiller():
    interface = solver_config._interface("copt")
    cobra.Configuration().solver = interface
    model = cobra.Model("incumbent_validation")
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    uptake = cobra.Reaction("UPTAKE", upper_bound=10)
    uptake.add_metabolites({a: 1})
    biomass = cobra.Reaction("BIOMASS")
    biomass.add_metabolites({b: -1})
    model.add_reactions([uptake, biomass])
    model.objective = biomass
    universal = cobra.Model("universal")
    conversion = cobra.Reaction("CONVERT")
    conversion.add_metabolites({a: -1, b: 1})
    universal.add_reactions([conversion])
    from cobra.flux_analysis.gapfilling import GapFiller

    gap = GapFiller(model, universal, demand_reactions=False, lower_bound=0.1)
    assert gap.model.optimize().status == "optimal"
    return gap


def test_gapfill_accepts_and_validates_resource_limited_incumbent(solved_copt_gapfiller):
    gap = solved_copt_gapfiller
    gap.model.solver._status = "time_limit"
    with patch.object(gap, "fill", side_effect=cobra.exceptions.OptimizationError(
            "Time limit")), patch.object(gap, "validate", wraps=gap.validate) as validate:
        reactions = fill_with_validated_incumbent(gap)

    assert [reaction.id for reaction in reactions] == ["CONVERT"]
    validate.assert_called_once_with(reactions)
    assert gap.original_model.slim_optimize() == 0
    gap.original_model.add_reactions(reactions)
    assert require_growth(gap.original_model, "Validated draft") == pytest.approx(10)


@pytest.mark.parametrize("status", ["infeasible", "nofeasible", "numeric", "in_progress"])
def test_gapfill_rejects_invalid_status_even_with_stale_primals(solved_copt_gapfiller, status):
    gap = solved_copt_gapfiller
    gap.model.solver._status = status
    with patch.object(gap, "fill", side_effect=cobra.exceptions.OptimizationError(
            "Infeasible")), patch.object(gap, "validate") as validate:
        with pytest.raises(RuntimeError, match="without a feasible incumbent"):
            fill_with_validated_incumbent(gap)
    validate.assert_not_called()


@pytest.mark.parametrize("error", [ValueError("unexpected adapter failure"),
                                 RuntimeError("unexpected adapter failure")])
def test_gapfill_preserves_unexpected_errors(solved_copt_gapfiller, error):
    gap = solved_copt_gapfiller
    with patch.object(gap, "fill", side_effect=error), patch.object(gap, "validate") as validate:
        with pytest.raises(type(error), match="unexpected adapter failure") as raised:
            fill_with_validated_incumbent(gap)
    assert raised.value is error
    validate.assert_not_called()


def test_gapfill_rejects_missing_incumbent_fluxes(solved_copt_gapfiller):
    gap = solved_copt_gapfiller
    gap.model.solver._status = "time_limit"
    with patch.object(gap, "fill", side_effect=cobra.exceptions.OptimizationError(
            "Time limit")), patch.object(type(gap.model.solver), "primal_values",
                                         new_callable=PropertyMock,
                                         return_value={"missing": math.nan}):
        with pytest.raises(RuntimeError, match="missing flux values"):
            fill_with_validated_incumbent(gap)


def test_gapfill_rejects_incumbent_without_validated_growth(solved_copt_gapfiller):
    gap = solved_copt_gapfiller
    gap.model.solver._status = "time_limit"
    with patch.object(gap, "fill", side_effect=cobra.exceptions.OptimizationError(
            "Time limit")), patch.object(gap.indicators[0], "_get_primal", return_value=0):
        with pytest.raises(RuntimeError, match="failed growth validation"):
            fill_with_validated_incumbent(gap)


@pytest.mark.parametrize("value", [None, math.nan, math.inf, -1, 0])
def test_growth_requires_finite_positive_objective(value):
    # Invalid objective values must be rejected before consulting native state.
    assert not has_growth(None, SimpleNamespace(objective_value=value))


@pytest.mark.parametrize("backend", ["gurobi", "copt"])
def test_gapfill_time_limit_does_not_repeat_native_solve(backend, monkeypatch):
    probe = subprocess.run([sys.executable, str(ROOT / "src/solver_config.py"),
                            "--probe", backend], capture_output=True, text=True, timeout=60)
    if probe.returncode:
        pytest.skip(f"{backend} unavailable: {probe.stderr.strip()}")
    interface = solver_config._interface(backend)
    solver = interface.Model("single_time_limited_solve")
    gap = SimpleNamespace(model=SimpleNamespace(solver=solver))
    configure_gapfiller(gap, 30, 8, 0.01)
    native = solver.problem
    get_param = native.getParam if backend == "copt" else lambda key: native.getParamInfo(key)[2]
    assert solver.configuration.presolve is True
    assert get_param("Presolve") == (1 if backend == "copt" else 2)

    if backend == "copt":
        from coptpy import COPT

        native_timeout_status = COPT.TIMEOUT
    else:
        import gurobipy

        native_timeout_status = gurobipy.GRB.TIME_LIMIT
    # Exercise the real optlang optimize and backend status mapping. Stub only
    # native execution so a deterministic resource limit never needs a long MIP.
    native_stub = SimpleNamespace(
        update=Mock(), reset=Mock(), solve=Mock(), optimize=Mock(),
        status=native_timeout_status,
        getAttr=Mock(return_value=native_timeout_status), setParam=Mock(),
        params=SimpleNamespace(Presolve=-1),
    )
    solver.update()
    monkeypatch.setattr(solver, "problem", native_stub)

    assert solver.optimize() == "time_limit"

    native_solve = native_stub.solve if backend == "copt" else native_stub.optimize
    native_solve.assert_called_once_with()
    assert solver.status == "time_limit"


@pytest.mark.parametrize("backend", ["gurobi", "copt"])
def test_real_gapfilling_and_reload(backend, tmp_path):
    # Probe independently so native license errors cannot kill pytest.
    probe = subprocess.run([sys.executable, str(ROOT / "src/solver_config.py"),
                            "--probe", backend], capture_output=True, text=True, timeout=60)
    if probe.returncode:
        pytest.skip(f"{backend} unavailable: {probe.stderr.strip()}")
    solver_config.select_solver(backend)
    model = cobra.Model("draft")
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    uptake = cobra.Reaction("UPTAKE")
    uptake.add_metabolites({a: 1})
    biomass = cobra.Reaction("BIOMASS")
    biomass.add_metabolites({b: -1})
    model.add_reactions([uptake, biomass])
    model.objective = biomass
    assert model.slim_optimize() == 0
    universal = cobra.Model("universal")
    conversion = cobra.Reaction("CONVERT")
    conversion.add_metabolites({a: -1, b: 1})
    conversion.gene_reaction_rule = "gene_a and gene_b"
    universal.add_reactions([conversion])
    from cobra.flux_analysis.gapfilling import GapFiller

    gap = GapFiller(model, universal, demand_reactions=False, lower_bound=0.1)
    configure_gapfiller(gap, 30, 8, 0.01)
    problem = gap.model.solver.problem
    get_param = problem.getParam if backend == "copt" else lambda key: problem.getParamInfo(key)[2]
    assert get_param("TimeLimit") == 30
    assert get_param("RelGap" if backend == "copt" else "MIPGap") == 0.01
    assert get_param("MemLimit" if backend == "copt" else "SoftMemLimit") == (8192 if backend == "copt" else 8)
    additions = fill_with_validation_retry(gap)
    assert [reaction.id for reaction in additions] == ["CONVERT"]
    model.add_reactions(additions)
    assert model.slim_optimize() > 0.1
    assert model.reactions.CONVERT.gene_reaction_rule == "gene_a and gene_b"
    optimum = model.slim_optimize()
    with model:
        cobra.manipulation.knock_out_model_genes(model, ["gene_a"])
        assert model.slim_optimize() == 0
    assert model.slim_optimize() == pytest.approx(optimum)
    parsimonious = cobra.flux_analysis.pfba(model)
    assert parsimonious.fluxes["BIOMASS"] == pytest.approx(optimum)
    path = tmp_path / "model.yml"
    cobra.io.save_yaml_model(model, str(path))
    reloaded = cobra.io.load_yaml_model(str(path))
    assert reloaded.solver.interface is model.solver.interface
    assert reloaded.slim_optimize() == pytest.approx(model.slim_optimize())
