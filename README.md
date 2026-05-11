# Modeling
The modeling was done on paper and the results can be found in the function sysID/EOMs/sysEOM .

Later experiments showed that the real plant had more effects to be modeled, hence: sysID/EOMs/sysEOM_asym .
Which accounts for Coulomb friction and also input deadzones and differential gains.

# Parameter estimation
Done using MATLAB's library for grey system identification.
A first experiment was conducted by keeping the cart steady: x = 0, but the full parameter estimation results in very close parameters and does not require too much time so the results of the "steady cart" experiment are unused (for reference estimation done in sysID a_fixedCartID.m).

The files **sysID/b_fullSysID.m** and **sysID/c_fullSys_validation.m** are for estimation and validation respectively.

# Linearization
Since the control techniques we have to use are linear, first we linearize around the unstable equilibrium.
The analytical linearization of the simple model is used for comparison with the more accurate one obtained using the simulink linearization tool on a model that uses the EOM in the funciton **sys_EOM_asym.m**.
The numerical solution has different values specially where the parameter 'c' appears in the analytical solution, this is probably because the dead-zones and Coulomb act "in parallel" to the cart friction to oppose motion. 