# Optimal Controllers for Load Restoration Problem to Improve Grid Resilience

Three Different Controllers: 
1.	Optimal Controller (Full Horizon Optimization) - Load  Restoration
2.	MPC - Load Restoration
3.	MPC - Python-Julia - Load Restoration

Notes:
1.	Clone the repo and put all the files in the repo in the same directory.
2.	The controllers have similar functions though their implementations are different.
3.	Please refer the DOC for the optimization problem definition and model formulation.
4.	Please change the path of the feeder data to your machine’s data path. The feeder OpenDSS file (IEEE13Nodeckt.dss) is in the folder “data”. 
5.	Please change the path of the exogenous dataset to your machine’s data path. The exogenous data (five_min_load_profile.csv) is in the folder “data”. 
6.	Controllers 1 and 2 are written in Julia and hence to run/use them it is required to install Julia/Ijulia and all the relevant Julia packages used.

      a.	Controller 1 is based on single-run (full-horizon) multi-step optimization

      b.	Controller 2 is based on multiple-run (for each step) single-step optimization, which is basically MPC except that, for the time being, it doesn’t look-ahead.

7.	Controller 3:

      a.	It is similar with Controller 2, but:

      b.	The data analytics (pre-processing and post-processing) is written in Python while the MPC (MPC.jl) is written in Julia and called from Python.

      c.	To run/use it, it is required to install Python/Ipython, Julia/Ijulia and all the relevant Python and Julia packages used and PyJulia

      d.	Typically, Julia is installed in /Applications, which isn't included in your PATH, and so the shell or conda-jupyter won't find it when you type Julia. That is PyJulia may not work with the default conda setting. To fix this: First find out the location of the Julia binary executable file in your machine. Then create a link to the executable and put it into the /usr/local/bin directory (which should already be in your path), so that typing julia is the exact equivalent of typing /Applications/Julia/.../julia. 
