# Model Predictive Control (MPC) for Load Restoration Problems in Active Power Distribution Systems

Four Different MPCs: 
1. MPC for copper plate system (without network flow and voltage constraints)
2. MPC for three-phase balanced system (with network flow and voltage constraints)
3. MPC for multi-phase unbalanced system (with network flow and voltage constraints)
4. MPC with reserve (co-optimization of power and reserve products)

Note:
The optimization models of the controllers are implemented using JuMP/Julia.
The data analytics (pre-processing and post-processing), renewable power forecasting and mpc run are implemented in Python.

To run the controllers:
1. Clone/download the repo
2. Install Python/Ipython, Julia/Ijulia, and all the relevant Python and Julia packages used.
3. Run the Python notebooks.
  
Note: 
Julia is installed in /Applications, which isn't included in your PATH, and so the shell or conda-jupyter won't find it when you type Julia. That is PyJulia may not work with the default conda setting. To fix this: First find out the location of the Julia binary executable file in your machine. Then create a link to the executable and put it into the /usr/local/bin directory (which should already be in your path), so that typing julia is the exact equivalent of typing /Applications/Julia/.../julia.