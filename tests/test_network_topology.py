"""
Distribution network topology tester.
"""

import pytest
import os
import warnings

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


class TestGridTopology:

    def test_grid_topology(self):

        bus = NETWORK[0]

        r = 0
        root_bus = None

        for b in bus.keys():
            
            l = len(bus[b].ancestor)
            
            if l > 1:
                warnings.warn('Network Not Radial; Bus ' f"{len(bus[b].index)}")
            elif l == 0:                
                root_bus = b
                r += 1   
                print("The root/substation bus is:", root_bus)                            
                
        if r == 0:
            warnings.warn("No root detected")
            root_bus = None
        elif r > 1:
            warnings.warn("More than one root detected - the network is not radial or it is multi-phase!")

        assert r == 1


# pytest test_network_topology.py