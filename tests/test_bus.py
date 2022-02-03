"""
Distribution grid bus tester.
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

class TestBus:

    def test_num_bus(self):

        bus = NETWORK[0]

        assert len(bus) == 13

    def test_slack_bus(self):

        bus = NETWORK[0]

        assert bus[1].is_root == True

    def test_simple_bus(self):

        bus = NETWORK[0]

        spec = np.array([10, False, 0.004911704881656805, 0.002600314349112426, 1.1025, 0.9025, 11, 9], dtype=object)

        assert bus[10].index == spec[0]
        assert  bus[10].is_root == spec[1]
        assert math.isclose(bus[10].d_P, spec[2], abs_tol=abs_tol)
        assert math.isclose(bus[10].d_Q, spec[3], abs_tol=abs_tol)  
        assert bus[10].v_max == spec[4]
        assert bus[10].v_min == spec[5]
        assert bus[10].children[0] == spec[6]
        assert bus[10].ancestor[0] == spec[7]


# pytest test_bus.py