# Overview
In this branch, vallina Lshaped method and stochastic decomposition (SD) are impelemented by MATLAB for two-stage stochastic LP programs, which read as

minimize    c'x + E[h(x, omega)]
subject to  A x (sense) b,
            x_lb <= x <= x_ub,
where "Ax (sense) b" can be mixed constraints consisting of "Ax = b", "Ax \leq b" as well as "Ax \geq b", and h(x, omega) = 
minimize    d‘y
subject to  D y (sense) r(omega) - C x,
            y_lb <= y <= y_ub.

**Requirements**
- Gruobi [https://www.gurobi.com]
- spUtilities/matlab_zyzhang [https://github.com/SMU-SODA/spUtilities/tree/matlab_zyzhang]

**Remarks**
- The uncertainty appears in the right-hand side r(omega) only.
- All scenario vectors are adjusted by subtracting the mean value, i.e., r(omega) - E[r(omega)].
- All possible scenario vectors and corresponding probabilities are generated explicitly in both Lshaped and SD. Please make sure that the problem size fits the memory.

**Ongoing Work**
- Randomness in the coefficient matrix and cost vector
- Efficient sampled paths management in SD for large-scale instances, e.g., ``ssn`` (over 1e70 scenarios).