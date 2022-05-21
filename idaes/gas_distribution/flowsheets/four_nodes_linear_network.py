import pyomo.environ as pyo
import idaes.core as idaes
import pyomo.network as network
#from pyomo.network.arc import Arc
from idaes.gas_distribution.properties.natural_gas import (
    NaturalGasParameterBlock,
)
from idaes.gas_distribution.unit_models.pipeline import GasPipeline
from idaes.gas_distribution.unit_models.compressor import (
    IsothermalCompressor as Compressor,
)
from idaes.gas_distribution.unit_models.node import PipelineNode
def four_nodes_linear():
    """
    The network looks something like this:
    s            s
    |            |
    0 (c)-> 1 -> 2 (c)-> 3
            |            |
            d            d
    
    """
    m = pyo.ConcreteModel()
    fs_config = {
        "dynamic": True,
        "time_set": [0.0, 20.0],
        "time_units": pyo.units.hr,
    }
    m.fs = idaes.FlowsheetBlock(default=fs_config)
    m.fs.properties = NaturalGasParameterBlock()

    pipeline_config = {
        "property_package": m.fs.properties,
        "finite_elements": 2,
    }
    m.fs.pipeline_set = pyo.Set(initialize=range(3))
    m.fs.pipeline = GasPipeline(m.fs.pipeline_set, default=pipeline_config)

    node_config = [
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 0,
            "n_outlet_pipelines": 1,
            "n_supplies": 1,
            "n_demands": 0,
        },
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 1,
            "n_supplies": 0,
            "n_demands": 1,
        },
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 1,
            "n_supplies": 1,
            "n_demands": 0,
        },
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 0,
            "n_supplies": 0,
            "n_demands": 1,
        },
    ]
    node_config = {i: config for i, config in enumerate(node_config)}
    m.fs.node_set = pyo.Set(initialize=range(4))
    m.fs.node = PipelineNode(m.fs.node_set, initialize=node_config)

    compressor_config = {"property_package": m.fs.properties}
    m.fs.compressor_set = pyo.Set(initialize=range(2))
    m.fs.compressor = Compressor(
        m.fs.compressor_set, default=compressor_config
    )

    # Connect compressors to pipelines
    # Should/could I make this easier?
    pipeline_idx_map = {0: 0, 1: 2}
    def compressor_to_pipeline_rule(fs, i):
        return (
            m.fs.compressor[i].outlet_port,
            m.fs.pipeline[pipeline_idx_map[i]].inlet_port,
        )
    m.fs.compressor_to_pipeline = network.Arc(
        m.fs.compressor_set, rule=compressor_to_pipeline_rule
    )

    # Note that we are adding a compressor, not a pipeline, to the
    # outlet here.
    m.fs.node[0].add_pipeline_to_outlet(m.fs.compressor[0])
    m.fs.node[1].add_pipeline_to_inlet(m.fs.pipeline[0])
    m.fs.node[1].add_pipeline_to_outlet(m.fs.pipeline[1])
    m.fs.node[2].add_pipeline_to_inlet(m.fs.pipeline[1])
    m.fs.node[2].add_pipeline_to_outlet(m.fs.compressor[1])
    m.fs.node[3].add_pipeline_to_inlet(m.fs.pipeline[2])

    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    expand_arcs.apply_to(m)


    ntfe = 2
    disc = pyo.TransformationFactory("dae.finite_difference")
    disc.apply_to(m, wrt=m.fs.time, nfe=ntfe)

    # Fix "design" variables
    for pipeline in m.fs.pipeline.values():
        pipeline.diameter.fix()
        pipeline.control_volume.length.fix()

    # Dynamic inputs:
    pred_dof = 9 * len(m.fs.time)
    # Initial conditions:
    pred_dof += (
        2 * (len(m.fs.pipeline[0].control_volume.length_domain) - 1)
    )
    pred_dof += (
        2 * (len(m.fs.pipeline[1].control_volume.length_domain) - 1)
    )
    pred_dof += (
        2 * (len(m.fs.pipeline[2].control_volume.length_domain) - 1)
    )

   

    # Fix predicted degrees of freedom:
    for t in m.fs.time:
        m.fs.node[0].supplies[0].state[t].mole_frac_comp[:].fix()
        m.fs.node[0].state[t].temperature.fix()
        m.fs.node[0].state[t].pressure.fix()
        m.fs.node[2].supplies[0].state[t].mole_frac_comp[:].fix()
        m.fs.node[2].supplies[0].state[t].flow_mol.fix()
        m.fs.node[1].demands[0].flow_mol[t].fix()
        m.fs.node[3].demands[0].flow_mol[t].fix()
        m.fs.compressor[0].boost_pressure[t].fix()
        m.fs.compressor[1].boost_pressure[t].fix()

    t0 = m.fs.time.first()
    x0 = m.fs.pipeline[0].control_volume.length_domain.first()
    xf = m.fs.pipeline[0].control_volume.length_domain.last()
    for x in m.fs.pipeline[0].control_volume.length_domain:
        # Here I assume that all three pipelines have the same
        # length domain.
        if x != x0:
            m.fs.pipeline[0].control_volume.pressure[t0, x].fix()
            m.fs.pipeline[1].control_volume.pressure[t0, x].fix()
            m.fs.pipeline[2].control_volume.pressure[t0, x].fix()
        if x != xf:
            m.fs.pipeline[0].control_volume.flow_mass[t0, x].fix()
            m.fs.pipeline[1].control_volume.flow_mass[t0, x].fix()
            m.fs.pipeline[2].control_volume.flow_mass[t0, x].fix()

    return m

def add_objective_to_model(
        m,
        dynamic=True,
        add_terminal_penalty=None,
        ):
    if add_terminal_penalty is None:
        add_terminal_penalty = dynamic
    m.supply_cost_coef = pyo.Param(initialize=0.0, mutable=True)
    m.electricity_cost_coef = pyo.Param(initialize=0.1, mutable=True)
    #m.demand_cost_coef = pyo.Param(initialize=1e6, mutable=True)
    # Demand cost coefficient needs to be scaled down for the units
    # of my model. The coefficient in the DAE code is 1e6 ((1e4 SCM)/hr)^-2.
    demand_conv_factor = (
        1e4
        * (0.72*pyo.units.kg / pyo.units.m**3)
        / (18.0*pyo.units.kg / pyo.units.kmol)
    )
    m.demand_cost_coef = pyo.Param(
        initialize=pyo.value(1e6*demand_conv_factor**-2),
        mutable=True,
    )
    #import pdb; pdb.set_trace()
    space = m.fs.pipeline[0].control_volume.length_domain
    x0 = space.first()
    xf = space.last()
    if dynamic:
        m.terminal_pressure_coef = pyo.Param(initialize=1e6, mutable=True)
        # Terminal flow penalty uses units of kg/hr instead of kmol/hr
        terminal_flow_conv_factor = 1e4 * (0.72*pyo.units.kg / pyo.units.m**3)
        m.terminal_flow_coef = pyo.Param(
            initialize=pyo.value(1e6*terminal_flow_conv_factor**-2),
            mutable=True,
        )
        m.terminal_pressure = pyo.Param(
            space,
            initialize=57.0,
            units=pyo.units.bar,
            mutable=True,
        )
        m.terminal_flow_mass = pyo.Param(
            space,
            initialize=16666.7,
            units=pyo.units.kmol/pyo.units.hr,
            mutable=True,
        )

    time = m.fs.time
    t0 = time.first()
    tf = time.last()
    supply = m.fs.node[0].supplies[0]
    demand = m.fs.node[1].demands[0]
    compressor = m.fs.compressor

    m.target_demand = pyo.Param(
        time,
        mutable=True,
        units=pyo.units.kmol/pyo.units.hr,
        initialize=16667.0,
    )

    if dynamic:
        m.supply_cost = pyo.Expression(expr=sum(
            m.supply_cost_coef * supply.flow_mol[t] * (t - time.prev(t))
            for t in time if t != t0
        ))

        m.boost_cost = pyo.Expression(expr=sum(
            m.electricity_cost_coef * compressor.power[t] * (t - time.prev(t))
            for t in time if t != t0
        ))

        m.demand_penalty = pyo.Expression(expr=sum(
            m.demand_cost_coef
            * (demand.flow_mol[t] - m.target_demand[t])**2
            * (t - time.prev(t))
            for t in time if t != t0
        ))

        m.terminal_pressure_cost = pyo.Expression(expr=sum(
            m.terminal_pressure_coef * (
                m.fs.pipeline.control_volume.pressure[tf, x]
                - m.terminal_pressure[x]
            )**2
            for x in space
        ))

        m.terminal_flow_mass_cost = pyo.Expression(expr=sum(
            m.terminal_flow_coef * (
                m.fs.pipeline.control_volume.flow_mass[tf, x]
                - m.terminal_flow_mass[x]
            )**2
            for x in space
        ))
    else:
        m.supply_cost = pyo.Expression(expr=sum(
            m.supply_cost_coef * supply.flow_mol[t] for t in time
        ))

        m.boost_cost = pyo.Expression(expr=sum(
            m.electricity_cost_coef * compressor.power[t] for t in time
        ))

        m.demand_penalty = pyo.Expression(expr=sum(
            m.demand_cost_coef
            * (demand.flow_mol[t] - m.target_demand[t])**2
            for t in time
        ))

    cost_expr = m.supply_cost + m.boost_cost + m.demand_penalty
    if dynamic and add_terminal_penalty:
        # TODO: Add terminal costs
        cost_expr += m.terminal_pressure_cost
        cost_expr += m.terminal_flow_mass_cost
    cost_expr *= 1e-6

    m.obj = pyo.Objective(expr=cost_expr)

m = four_nodes_linear()
pyo.display(m)
add_objective_to_model(m)