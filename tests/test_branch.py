"""
Distribution grid branch/line tester.
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


class TestBranch:

    def test_num_branch(self):

        branch = NETWORK[1]

        assert len(branch) == 12
    
    def test_simple_branch(self):

        branch = NETWORK[1]

        spec = np.array([6, 12, 8, 0.114, 0.179, 3.974509847458756], dtype=object)

        assert branch[6].index == spec[0]
        assert branch[6].to_node == spec[1]
        assert branch[6].from_node == spec[2]
        assert math.isclose(branch[6].r, spec[3], abs_tol=abs_tol)
        assert math.isclose(branch[6].x, spec[4], abs_tol=abs_tol)
        assert math.isclose(branch[6].b, spec[5], abs_tol=abs_tol)
      

# pytest test_branch.py