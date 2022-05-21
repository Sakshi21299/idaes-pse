# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 14:14:22 2022

@author: ssnaik
"""

#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
"""
import pyomo.environ as pyo
from collections import namedtuple
from pyomo.dae.flatten import flatten_dae_components
from pyomo.util.subsystems import TemporarySubsystemManager

from idaes.gas_distribution.flowsheets.simple_network_model import (
    make_simple_model,
    add_objective_to_model,
)

from idaes.apps.nmpc.dynamic_data import (
    interval_data_from_time_series,
    load_inputs_into_model,
)
from idaes.core.util.model_statistics import variables_near_bounds_generator
import matplotlib.pyplot as plt


def _plot_time_indexed_data(
        data,
        names,
        show=True,
        prefix=None,
        ):
    """ 
    data:
    (
        [t0, ...],
        {
            str(cuid): [value0, ...],
        },
    )
    names: list of str(cuids) that we will actually plot.
    """
    if prefix is None:
        prefix = ""
    plt.rcParams.update({"font.size": 16})
    time, name_map = data
    for i, name in enumerate(names):
        values = name_map[name]
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(time, values, linewidth=3)

        # Logic for adding useful names and axis labels to plots we care
        # about using for other purposes.
        ax.set_xlabel("Time (hr)")
        ax.set_title(name)
        fig.tight_layout()
        if show:
            plt.show()
        else:
            fig.savefig(prefix + "state%s.png" % i, transparent=True)


"""
A data structure we'll use for time series data from Pyomo variables.
Just a list of time points, followed by a dict mapping ComponentUIDs to
lists of variable values. I.e.:
(
    [t0, t1, ...],
    {
        ComponentUID(var): [val0, val1, ...],
    },
)
We implicitly assume that each variable/CUID has the same number of
values as there are time points.
"""
TimeSeriesTuple = namedtuple("TimeSeriesTuple", ["time", "data"])


def get_nominal_demand_sequence():
    """
    Get the nominal piecewise constant sequence of demand values.
    """
    sample_points = [0.0, 4.0, 20.0]
    # This is the "nominal demand trajectory"
    # I don't know what I was thinking by using different units here.
    # TODO: Change this later and actually use unit-ed expressions here.
    # 30 kg/hr -> kmol/hr
    val0 = pyo.value(30.0 * 1e4 / 18.0)
    # 50 SCM/hr -> kmol/hr
    val1 = pyo.value(50.0 * 1e4 * 0.72 / 18.0)
    demand_cuid = pyo.ComponentUID("target_demand")
    # This data is meant to be a piecewise constant demand trajectory.
    # val0 between the first two sample points, and val1 between the last
    # two sample points.
    # At a boundary, the value of the interval to the right is used.
    demand_data = TimeSeriesTuple(
        sample_points,
        {
            demand_cuid: [val0, val0, val1],
        },
    )
    return demand_data


def get_perturbed_demand_sequence():
    """
    Get the piecewise constant sequence of perturbed demand values.
    """
    sample_points = [0.0, 4.0, 20.0]
    val0 = 30.0 * 1e4 / 18.0
    val1 = 50.0001 * 1e4 * 0.72 / 18.0
    demand_cuid = pyo.ComponentUID("target_demand")
    demand_data = TimeSeriesTuple(
        sample_points,
        {
            demand_cuid: [val0, val0, val1],
        },
    )
    return demand_data


def run_dynamic_optimization():
    """
    Just run a dynamic optimization problem and return variable data
    from the solve.
    """
    demand = get_nominal_demand_sequence()
    m = solve_boost_pressure_given_demand(demand)
    time = m.fs.time

    #
    # Extract values we care about and plot
    #
    scalar_vars, dae_vars = flatten_dae_components(m, time, pyo.Var)
    sim_data = (
        list(time),
        {
            str(pyo.ComponentUID(var.referent)): [var[t].value for t in time]
            for var in dae_vars
        },
    )
    return sim_data


def solve_boost_pressure_given_demand(demand):
    """
    Create and solve a model in which boost pressures are degrees of freedom
    used to satisfy the given demand data. The optimization problem also
    includes a terminal penalty penalizing deviation from a steady state with
    demand at the end of the horizon, and a lower bound on pressure at the
    demand node.

    Parameters
    ----------
    demand: TimeSeriesTuple
        A tuple where the first entry is a list of time points and the
        second entry is a dict mapping CUIDs to lists of values. The
        data correspond to a piecewise-constant sequence of inputs.

    """
    nxfe = 4

    ipopt = pyo.SolverFactory("ipopt")

    demand_cuid = pyo.ComponentUID("target_demand")

    #
    # Create steady state model for initial conditions and extract data
    #
    m = make_simple_model(nxfe=nxfe, dynamic=False)
    prop = m.fs.properties
    mw = prop.natural_gas.mw
    # Inlet flow of the initial condition steady state is the same as
    # the initial target demand. This is because at steady state the flow
    # rate is constant along the pipeline.
    inlet_flow = demand.data[demand_cuid][0]
    t0 = m.fs.time.first()
    m.fs.nodes[0].supplies[0].flow_mol[t0].fix(inlet_flow)
    m.fs.nodes[0].state[t0].pressure.fix(50.0*pyo.units.bar)
    m.fs.compressor.boost_pressure[t0].fix(7.0*pyo.units.bar)
    res = ipopt.solve(m, tee=True)
    pyo.assert_optimal_termination(res)
    scalar_vars, dae_vars = flatten_dae_components(m, m.fs.time, pyo.Var)
    initial_scalar_data = {
        str(pyo.ComponentUID(var)): var.value
        for var in scalar_vars
    }
    initial_data = {
        str(pyo.ComponentUID(var.referent)): var[t0].value
        for var in dae_vars
    }
    ###

    p_demand_min = 50.5
    #
    
    #
    add_objective_to_model(m, dynamic=False)
    m.fs.nodes[0].supplies[0].flow_mol[t0].unfix()
    m.fs.compressor.boost_pressure[t0].unfix()
    m.fs.compressor.boost_pressure[t0].setlb(0.0)
    m.fs.compressor.boost_pressure[t0].setub(100.0)
    m.fs.nodes[1].state[t0].pressure.setlb(50.5)
    m.fs.nodes[1].demands[0].flow_mol[t0].unfix()

    # Target demand for the terminal steady state is demand at the end
    # of the provided sequence.
    target_demand = demand.data[demand_cuid][-1]
    m.target_demand[:].set_value(target_demand)
    res = ipopt.solve(m, tee=True)
    pyo.assert_optimal_termination(res)
    target_data = {
        str(pyo.ComponentUID(var.referent)): var[t0].value
        for var in dae_vars
    }
    ###

    #
    # Construct dynamic model
    #
    horizon = 20.0
    ntfe = 20
    time_scheme = "BACKWARD"
    m = make_simple_model(
        horizon=horizon,
        ntfe=ntfe,
        nxfe=nxfe,
        dynamic=True,
        time_scheme=time_scheme,
    )
    add_objective_to_model(m, dynamic=True)
    time = m.fs.time
    t0 = time.first()
    ###

    #
    # Set terminal costs from target steady state
    #
    space = m.fs.pipeline.control_volume.length_domain
    for x in space:
        m.terminal_pressure[x].set_value(
            # Ideally I would have a map to the variable in the dynamic model,
            # then a map from that variable to its terminal parameter.
            target_data["fs.pipeline.control_volume.pressure[*,%s]" % x]
        )
        m.terminal_flow_mass[x].set_value(
            target_data["fs.pipeline.control_volume.flow_mass[*,%s]" % x]
        )

    #
    # Fix degrees of freedom
    #
    m.fs.compressor.boost_pressure[:].fix()
    m.fs.nodes[0].state[:].pressure.fix()
    m.fs.nodes[1].demands[0].flow_mol[:].fix()
    ###

    #
    # Initialize dynamic model with steady state data
    #
    for name, val in initial_data.items():
        var = m.find_component(name)
        for t in m.fs.time:
            var[t].set_value(val)
    for name, val in initial_scalar_data.items():
        var = m.find_component(name)
        var.set_value(val)
    m.fs.pipeline.control_volume.material_accumulation[...].set_value(0.0)
    m.fs.pipeline.control_volume.flow_mass_dt[:,:].set_value(0.0)
    ###

    #
    # Sanity check solve - should have zero degrees of freedom and no
    # infeasibility
    #
    res = ipopt.solve(m, tee=True)
    pyo.assert_optimal_termination(res)

    #
    # Construct target demand sequence and load into model
    #
    # Note that here we are setting the value of a mutable parameter
    input_interval_data = interval_data_from_time_series(demand)
    load_inputs_into_model(m, time, input_interval_data)
    ###

    #
    # Make sure proper degrees of freedom are unfixed
    #
    m.fs.nodes[0].state[:].pressure.fix()
    m.fs.compressor.boost_pressure[:].setub(100.0)
    for t in time:
        if t != t0:
            m.fs.compressor.boost_pressure[t].unfix()
            m.fs.nodes[1].demands[0].flow_mol[t].unfix()

            m.fs.nodes[1].state[t].pressure.setlb(p_demand_min)

    m.fs.pipeline.control_volume.area.fix()
    m.fs.pipeline.diameter_eqn.deactivate()
    res = ipopt.solve(m, tee=True)
    pyo.assert_optimal_termination(res)
    var_near_bounds = variables_near_bounds_generator(m, tol = 10**-4)
    var_below_bound = []
    for var in var_near_bounds:
        if pyo.value(var)< 50.5:
            print(var)
            var_below_bound.append(var)
    return m


def solve_given_boost_pressure_with_new_demand(m, demand):
    """
    Solve a dynamic optimization problem with boost pressure fixed to their
    current values, using the demand specified. Terminal penalties are
    eliminated by setting their coefficients to zero. These changes are
    only made temporarily so that we don't have to remember anything about
    the "state" of our model (other than variable values) in calling functions.

    Parameters
    ----------
    m: ConcreteModel
    demand: TimerSeriesTuple

    """
    to_reset = (
        list(m.terminal_pressure_coef.values())
        + list(m.terminal_flow_coef.values())
        + list(m.target_demand.values())
    )
    to_fix = list(m.fs.compressor.boost_pressure.values())

    # Just do this here so I can initialize to initial conditions
    time = m.fs.time
    t0 = time.first()
    scalar_vars, dae_vars = flatten_dae_components(m, time, pyo.Var)

    with TemporarySubsystemManager(
            to_fix=to_fix, to_reset=to_reset
            ):
        m.terminal_pressure_coef[:].set_value(0.0)
        m.terminal_flow_coef[:].set_value(0.0)

        # Initialize to initial conditions.
        # The previous solve (whatever that is) may be better initialization,
        # but I do this for consistency.
        for var in dae_vars:
            for t in time:
                # We do not want to override boost pressures, which are fixed.
                if not var[t].fixed:
                    var[t].set_value(var[t0].value)

        # Load new demand. These will get overridden by the original
        # values after the solve. Not sure how I feel about this.
        input_interval_data = interval_data_from_time_series(demand)
        load_inputs_into_model(m, time, input_interval_data)

        # Have some infeasibility here as I have initialized to initial
        # conditions.
        ipopt = pyo.SolverFactory("ipopt")
        res = ipopt.solve(m, tee=True)
        pyo.assert_optimal_termination(res)
        import pdb; pdb.set_trace()
        var_near_bounds = variables_near_bounds_generator(m, tol = 10**-4)
        var_below_bound = []
        for var in var_near_bounds:
            print(pyo.value(var))
            if pyo.value(var)< 50.5:
                print(var)
                var_below_bound.append(var)
        #The pressure at node[1] is exactly at it's bound after 5 hours. 
        
        


def main():
    demand = get_nominal_demand_sequence()
    m = solve_boost_pressure_given_demand(demand)
    time = m.fs.time

    # Flatten variables
    scalar_vars, dae_vars = flatten_dae_components(m, time, pyo.Var)
    # Hard-code names of variables we care about
    names_to_plot = [
        "fs.nodes[0].state[*].flow_mol",
        "fs.nodes[1].state[*].flow_mol",
        "fs.compressor.boost_pressure[*]",
        "fs.nodes[1].state[*].pressure",
    ]

    # Extract data from nominal optimization problem
    sim_data = (
        list(time),
        {
            str(pyo.ComponentUID(var.referent)): [var[t].value for t in time]
            for var in dae_vars
        },
    )
    _plot_time_indexed_data(
        sim_data, names_to_plot, prefix="nominal_", show=False
    )

    # Extract data from "optimization problem" with perturbed demand
    new_demand = get_perturbed_demand_sequence()
    solve_given_boost_pressure_with_new_demand(m, new_demand)
    sim_data = (
        list(time),
        {
            str(pyo.ComponentUID(var.referent)): [var[t].value for t in time]
            for var in dae_vars
        },
    )
    _plot_time_indexed_data(
        sim_data, names_to_plot, prefix="perturbed_", show=False
    )


if __name__ == "__main__":
    main()
