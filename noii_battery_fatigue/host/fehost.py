"""The benchmark host's element interface, on the public abaqus_ufl.fe (github.com/tengzhang48/abaqus_ufl).

Same functions as the host the agents were given (load_contract, build_module, ContractUEL,
residual_weights, solve). Building, the f2py wrapper, the element adapter and the solver are public
abaqus_ufl.fe: build_external_uel, CompiledAbaqusElement.from_module, CompiledElementState and
solve_steps_staggered. Only the benchmark contract stays here: parsing uel_contract.json, the DOF
names, and the task's promised call convention TIME(1) = TIME(2) = END of the increment and
KSTEP = KINC = 1.
"""
import json
import os

from abaqus_ufl.fe import (CompiledAbaqusElement, CompiledElementState, IncrementTooLarge,  # noqa: F401
                           build_external_uel, solve_steps_staggered)

REQUIRED_DOFS = ("ux", "uy", "c", "d")


def load_contract(workspace):
    path = os.path.join(workspace, "element", "uel_contract.json")
    with open(path) as f:
        c = json.load(f)
    names = c["dof_names"]
    if len(names) != int(c["dof_per_node"]):
        raise ValueError("dof_names must list exactly dof_per_node names")
    for n in REQUIRED_DOFS:
        if n not in names:
            raise ValueError(f"dof_names must include {REQUIRED_DOFS}; missing {n!r}")
    return c


def build_module(workspace, contract, workdir, module_name):
    """Compile the contract's UEL sources (in order) and link them with the public wrapper."""
    sources = contract["source"] if isinstance(contract["source"], list) else [contract["source"]]
    return build_external_uel([os.path.join(os.path.abspath(workspace), s) for s in sources],
                              workdir=workdir, module_name=module_name,
                              compile_flags=list(contract.get("compile_flags", [])))


class EndOfIncrementUEL(CompiledAbaqusElement):
    """The task's call contract (README section 1): TIME(1) = TIME(2) = END of the increment, KINC = 1."""

    def begin_increment(self, kinc, time_begin=None, dt=None):
        super().begin_increment(1, time_begin=time_begin, dt=dt)
        self.time = self.time + self.dt

    def comp(self, name):
        return self.dof_names.index(name)


def ContractUEL(module, contract, n_elem, dt):
    """A compiled Quad4 UEL set up from uel_contract.json (run_problem.py calls this)."""
    element = EndOfIncrementUEL.from_module(
        module, dof_per_node=int(contract["dof_per_node"]), props=contract["props"],
        n_elem=n_elem, nsvars=int(contract.get("nsvars", 0)),
        svars_initial=contract.get("svars_initial"), jprops=contract.get("jprops", []),
        lflags=contract.get("lflags", [1, 0, 1, 0, 0, 0]), dt=dt)
    element.dof_names = list(contract["dof_names"])
    return element


def residual_weights(element, contract, n_nodes=None):
    return {element.comp(name): 1.0 / float(scale)
            for name, scale in (contract.get("residual_scale") or {}).items()}


def solve(nodes, elems, element, dirichlet_fn, n_steps, dt, *, U0=None, tol=1e-8,
          max_cutbacks=10, residual_scale=None, step_callback=None, commit_fn=None,
          batch_fn=None, max_field_newton=25, verbose=False):
    """Monolithic transactional solve of all DOFs as one field, with cutbacks."""
    dpn = element.dof_per_node
    return solve_steps_staggered(
        nodes, elems, dpn, element.element_rk, dirichlet_fn, [list(range(dpn))], n_steps,
        batch_fn=batch_fn or element.element_rk_batch,
        commit_fn=commit_fn or element.commit_step,
        state=CompiledElementState([element]), dt=dt, tol=tol, max_cutbacks=max_cutbacks,
        cutback_on_element_rejection=True, residual_scale=residual_scale, U0=U0,
        step_callback=step_callback, max_field_newton=max_field_newton, verbose=verbose,
        return_report=True)
