"""Legacy fallback for optlang-copt 0.2.0; 0.2.1 fixes these issues natively."""


def apply_copt_compatibility(interface):
    """Keep mass-balance expressions consistent when temporary variables leave.

    COBRA sets new reaction coefficients directly on existing constraints.
    optlang-copt 0.2.0 snapshots those native rows on variable deletion, but
    optlang's variable-to-constraint index may not include the edited rows.
    Linear models remove native variables directly; quadratic models retain the
    adapter's defensive rebuild, with all removed symbols substituted first.
    This uses the installed adapter; no site-packages files are changed.
    """
    from importlib.metadata import version
    from math import nan
    from optlang.interface import Model

    if version("optlang-copt") != "0.2.0":
        return
    if getattr(interface.Model, "_alphagem_rollback_compatibility", False):
        return

    set_native_objective = interface.Model._set_native_objective
    get_primal_values = interface.Model._get_primal_values
    get_objective_value = interface.Objective.value.fget
    get_model_state = interface.Model.__getstate__
    set_model_state = interface.Model.__setstate__

    def model_state_with_configuration(self):
        state = get_model_state(self)
        state['_alphagem_configuration'] = self.configuration.__getstate__()
        return state

    def restore_model_configuration(self, state):
        set_model_state(self, state)
        configuration = state.get('_alphagem_configuration')
        if configuration is not None:
            self.configuration.__setstate__(configuration)

    def objective_value_with_missing(self):
        if self.problem is None:
            return None
        native = self.problem.problem
        has_solution = native.hasmipsol if self.problem.is_integer else native.haslpsol
        if (not has_solution or self.problem.status in
                {"infeasible", "unbounded", "infeasible_or_unbounded"}):
            # COPT's 1e30 no-solution sentinel must never appear as growth.
            return nan
        return get_objective_value(self)

    def primal_values_with_missing(self):
        # COBRA subtracts forward/reverse primals even for infeasible drafts.
        # Missing values must remain numeric missing data, never feasible fluxes.
        return [nan if value is None else value for value in get_primal_values(self)]

    def track_quadratic_objective(self, expression, direction):
        if self.problem.HasQObj:
            self._alphagem_had_quadratic_objective = True
            # The adapter's setObjective does not replace existing QP terms.
            self.problem.delQuadObj()
        set_native_objective(self, expression, direction)
        if self.problem.HasQObj:
            self._alphagem_had_quadratic_objective = True

    def remove_variables(self, variables):
        self.objective._expression = self.objective.expression
        self.objective._expression_expired = False
        if not (self.problem.HasQObj or
                getattr(self, "_alphagem_had_quadratic_objective", False)):
            # Gene deletion can remove thousands of reaction variables. Building
            # a fresh native model for each deletion makes draft transfer unusable.
            native_variables = [variable._internal_variable for variable in variables]
            Model._remove_variables(self, variables)
            self.problem.remove(native_variables)
            self.problem.update()
            return
        for constraint in self._constraints:
            constraint._expression = constraint.expression
        Model._remove_variables(self, variables)
        replacements = {variable: 0 for variable in variables}
        for constraint in self._constraints:
            constraint._expression = constraint._expression.xreplace(replacements)
        self._rebuild_native_problem()

    interface.Model._set_native_objective = track_quadratic_objective
    interface.Model._get_primal_values = primal_values_with_missing
    interface.Model.__getstate__ = model_state_with_configuration
    interface.Model.__setstate__ = restore_model_configuration
    interface.Objective.value = property(objective_value_with_missing)
    interface.Model._remove_variables = remove_variables
    interface.Model._alphagem_rollback_compatibility = True
