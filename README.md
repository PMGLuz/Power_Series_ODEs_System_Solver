# Singular_ODEs_System_Solver
This mathematica package solves a family of systems of linear ODEs with an isolated singular point of the first kind, finding solutions in a power series form around the singular point.


Given a system of ordinary differential equations in a single variable of the form:
```
$$
\frac{d \mathsf{W}}{dx} = \left(\frac{1}{x} \mathrm{R} + \Theta \right) \mathsf{W}\,,
$$
```
where $\mathrm{R}$ is a constant matrix and $\Theta$ is an analytic matrix, this package finds the intermediate matrices for the solutions using a power series expansion, around the singular point.
The method can be found in chapter 4 of E. A. Coddington and N. Levinson, *Theory of ordinary differential equations* (McGraw-Hill, New York, 1955) or chapter 6.2 of E. A. Coddington an R. Carlson, *Linear Ordinary Differential Equations* (Society for Industrial and Applied Mathematics, Philadelphia, 1997).

A mathematical notebook with examples is provided. Some of the provided examples show how the method can be used to find the analytical eigenfunctions of the radial perturbations of static, spherically symmetric spacetimes in the theory of General Relativity.   
