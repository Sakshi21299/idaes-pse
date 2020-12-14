##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Python script to read costing components

"""
__author__ = "Costing Team (A. Noring and M. Zamarripa)"
__version__ = "1.0.0"

import os
import json

full_path = os.path.realpath(__file__)
directory = os.path.dirname(full_path)

''' The costing exponents dictionary contains information from the QGESS on
capital cost scaling methodology (DOE/NETL-2019/1784). Specifically it includes
scaling exponents, valid ranges for the scaled parameter, and units for those
ranges. It is important to note the units only apply to the ranges and are not
neccessarily the units that the reference parameter value will be given in.
This dictionary is nested with the following structure:
    tech type --> account --> property name --> property value'''

with open(directory + "\BB_costing_exponents.json", 'r') as file:
    BB_costing_exponents = json.load(file)


'''
The costing params dictionary contains information from the BBR4 COE
spreadsheet. It includes the total plant cost (TPC), reference parameter value,
and units for that value.

Some accounts are costed using two different reference
parameters, these accounts have been divided into two separate accounts
following
the naming convention x.x.a and x.x.b.
This dictionary is nested with the following structure:
tech type --> CCS --> account --> property name --> property values'''
with open(directory + "\BB_costing_parameters.json", 'r') as file:
    BB_costing_params = json.load(file)


with open(directory + "\sCO2_costing_parameters.json", 'r') as file:
    sCO2_costing_params = json.load(file)
