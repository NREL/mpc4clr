"""
Distributed energy resource (DER) tester.
"""

import pytest
import numpy as np
import os
import math

import julia
from julia.api import Julia

jl = Julia(compiled_modules=False)  
jl = julia.Julia()

current_dir = os.getcwd()
os.chdir(os.path.dirname(current_dir) + "/mpc_three_phase_balanced_system")
    
jl.include("data_handler_threephase.jl")
jl.include("lr_threephase_opf_mpc.jl")

NETWORK = jl.load_case_data(datafile="13buscase")

abs_tol = 1e-05


class TestDER:

    def test_wind(self):

        wind = NETWORK[3]

        assert len(wind) == 1

        spec = np.array([11, 0.023113905325443787, 0.01733542899408284, 0.028892381656804734], dtype=object)
      
        assert wind["wt1"].bus_idx == spec[0]
        assert math.isclose(wind["wt1"].w_P_max, spec[1], abs_tol=abs_tol)
        assert math.isclose(wind["wt1"].w_Q_max, spec[2], abs_tol=abs_tol)
        assert math.isclose(wind["wt1"].w_S_max, spec[3], abs_tol=abs_tol)

    def test_solar(self):
        
        pv = NETWORK[4]

        assert len(pv) == 1

        spec = np.array([13, 0.01733542899408284, 0.013001571745562131, 0.02166928624260355], dtype=object)

        assert pv["pv1"].bus_idx == spec[0]
        assert math.isclose(pv["pv1"].p_P_max, spec[1], abs_tol=abs_tol)
        assert math.isclose(pv["pv1"].p_Q_max, spec[2], abs_tol=abs_tol)
        assert math.isclose(pv["pv1"].p_S_max, spec[3], abs_tol=abs_tol)

    def test_microturbine(self):

        microturbine = NETWORK[2]

        assert len(microturbine) == 1

        spec = np.array([1, 0.023113905325443787, 0.01733542899408284], dtype=object)

        assert microturbine["g1"].bus_idx == spec[0]        
        assert math.isclose(microturbine["g1"].g_P_max, spec[1], abs_tol=abs_tol)
        assert math.isclose(microturbine["g1"].g_Q_max, spec[2], abs_tol=abs_tol)

    def test_battery(self):

        battery = NETWORK[5]

        assert len(battery) == 1

        spec = np.array([9, 0.011556952662721894, 0.00866771449704142, 100.0, 20.0, 1.0], dtype=object)

        assert battery['s1'].bus_idx == spec[0]        
        assert math.isclose(battery["s1"].s_P_max, spec[1], abs_tol=abs_tol)
        assert math.isclose(battery["s1"].s_Q_max, spec[2], abs_tol=abs_tol)
        assert battery["s1"].s_SOC_max == spec[3]  
        assert battery["s1"].s_SOC_min == spec[4]         
        assert battery["s1"].s_cap == spec[5]
        
          
# pytest test_DER.py