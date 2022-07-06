# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 11:22:34 2022

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
This file contains utilities for constructing a flowsheet for Kai's 5 demands,
3 compressors, 1 supply gas pipeline network.
"""
import pyomo.environ as pyo
import idaes.core as idaes
import pyomo.network as network
from idaes.models_extra.gas_distribution.properties.natural_gas import (
    NaturalGasParameterBlock,
)
from idaes.models_extra.gas_distribution.unit_models.pipeline import GasPipeline
from idaes.models_extra.gas_distribution.unit_models.compressor import (
    IsothermalCompressor as Compressor,
)
from idaes.models_extra.gas_distribution.unit_models.node import PipelineNode

from idaes.core.util.model_statistics import degrees_of_freedom

def make_model(
    dynamic=True,
    scenario = 1,
    nxfe=2,
    space_method="dae.finite_difference",
    space_scheme="FORWARD",
    ntfe=40,
    horizon=20.0,
    time_method="dae.finite_difference",
    time_scheme="BACKWARD",
    ):
    
    m = pyo.ConcreteModel()
    default = {"dynamic": dynamic}
    if dynamic:
        default["time_set"] = [0.0, horizon]
        default["time_units"] = pyo.units.hr
    
    m.fs = idaes.FlowsheetBlock(default=default)
    m.fs.properties = NaturalGasParameterBlock()


    

    node_configs = [
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
            "n_outlet_pipelines": 2,
            "n_supplies": 0,
            "n_demands": 0,
        },
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 2,
            "n_supplies": 0,
            "n_demands": 1,
        },
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 0,
            "n_supplies": 0,
            "n_demands": 1,
        },
        
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 0,
            "n_supplies": 0,
            "n_demands": 1,
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
            "n_outlet_pipelines": 0,
            "n_supplies": 0,
            "n_demands": 1,
        }
    ]
    m.fs.node_set = pyo.Set(initialize=list(range(len(node_configs))))
    node_configs = {i: config for i, config in enumerate(node_configs)}
    m.fs.nodes = PipelineNode(m.fs.node_set, initialize=node_configs)
   
    pipeline_config = {
        "property_package": m.fs.properties,
        "finite_elements": nxfe,
        "transformation_method": space_method,
        "transformation_scheme": space_scheme,
        "has_holdup": True,
    }
    m.fs.pipeline_set = pyo.Set(initialize=range(6))
    m.fs.pipeline = GasPipeline(m.fs.pipeline_set, default=pipeline_config)
    
    compressor_config = {"property_package": m.fs.properties}
    m.fs.compressor_set = pyo.Set(initialize=range(3))
    m.fs.compressor = Compressor(
        m.fs.compressor_set, default=compressor_config
    )
    m.fs.compressor[:].scenario_multiplication_parameter[:] = scenario
    

    # Connect compressors to pipelines

    m._compressor_to_pipeline_1 = network.Arc(
        ports=(m.fs.compressor[0].outlet_port,
        m.fs.pipeline[0].inlet_port)
    )
    
    m._compressor_to_pipeline_2 = network.Arc(
        ports=(m.fs.compressor[1].outlet_port,
        m.fs.pipeline[1].inlet_port)
    )
    
    m._compressor_to_pipeline_3 = network.Arc(
        ports=(m.fs.compressor[2].outlet_port,
        m.fs.pipeline[4].inlet_port)
    )
    
    m.fs.nodes[0].add_pipeline_to_outlet(m.fs.compressor[0])
    
    m.fs.nodes[1].add_pipeline_to_inlet(m.fs.pipeline[0])
    m.fs.nodes[1].add_pipeline_to_outlet(m.fs.compressor[1])
    m.fs.nodes[1].add_pipeline_to_outlet(m.fs.compressor[2])
    
    m.fs.nodes[2].add_pipeline_to_inlet(m.fs.pipeline[1])
    m.fs.nodes[2].add_pipeline_to_outlet(m.fs.pipeline[2])
    m.fs.nodes[2].add_pipeline_to_outlet(m.fs.pipeline[3])
    
    m.fs.nodes[3].add_pipeline_to_inlet(m.fs.pipeline[2])
    m.fs.nodes[4].add_pipeline_to_inlet(m.fs.pipeline[3])
    
    m.fs.nodes[5].add_pipeline_to_inlet(m.fs.pipeline[4])
    m.fs.nodes[5].add_pipeline_to_outlet(m.fs.pipeline[5])
    
    m.fs.nodes[6].add_pipeline_to_inlet(m.fs.pipeline[5])
    
    
    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    expand_arcs.apply_to(m)
    
    #Set pipeline length
    cv = m.fs.pipeline[:].control_volume
    m.fs.pipeline[:].diameter.fix(0.92*pyo.units.m)
    cv.length.fix(300*pyo.units.km)
    
    
    # Initial conditions:
   
    x0 = m.fs.pipeline[0].control_volume.length_domain.first()
    xf = m.fs.pipeline[0].control_volume.length_domain.last()
    j = next(iter(m.fs.properties.component_list))
    if dynamic:
        t0 = m.fs.time.first()
        for x in m.fs.pipeline[0].control_volume.length_domain:
            # Here I assume that all three pipelines have the same
            # length domain.
            if x != x0:
                m.fs.pipeline[0].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[1].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[2].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[3].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[4].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[5].control_volume.pressure[t0, x].fix()
            if x != xf:
                m.fs.pipeline[0].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[1].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[2].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[3].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[4].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[5].control_volume.flow_mass[t0, x].fix()
            
        cv.momentum_balance[t0, xf].deactivate()
        
        disc = pyo.TransformationFactory(time_method)
        disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme=time_scheme)

            
    # Fix "dynamic inputs." This needs to be done after a potential
    # discretization transformation.

    m.fs.nodes[0].state[:].mole_frac_comp[j].fix()
    m.fs.nodes[0].state[:].temperature.fix()
  
    
    return m 
    
def make_steady_model(
        scenario=1,
        nxfe=4,
        to_fix=None,
        input_data=None,
        tee=True,
        ):
    if to_fix is None:
        to_fix = []
    if input_data is None:
        input_data = {}
    m = make_model(scenario = scenario, nxfe=nxfe, dynamic=False)

    for cuid in to_fix:
        var = m.find_component(cuid)
        var[:].fix()

    for cuid, val in input_data.items():
        var = m.find_component(cuid)
        var[:].set_value(val)
    
    # Ipopt bound multipliers (obtained from solution)
    # m.ipopt_zL_out = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    # m.ipopt_zU_out = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    # # Ipopt bound multipliers (sent to solver)
    # m.ipopt_zL_in = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    # m.ipopt_zU_in = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    # m.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT_EXPORT)
    
    
    
    ipopt = pyo.SolverFactory("ipopt")
    res = ipopt.solve(m, tee=tee)
    pyo.assert_optimal_termination(res)

    return m

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
