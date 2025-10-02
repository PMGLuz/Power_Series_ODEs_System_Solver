# Power_Series_ODEs_System_Solver
The Mathematica package solves a family of systems of linear ODEs with an isolated singular point of the first kind, finding solutions in a power series form around the singular point.


Let a system of $n$ ordinary differential equations for n unknowns of the form:
```math
\frac{d \mathsf{W}}{dx} = \left(\frac{1}{x} \mathrm{R} + \Theta \right) \mathsf{W}\,,
```
where $\mathrm{R}$ is a constant $n\times n$ real matrix, and $\Theta$ is a $n\times n$ real analytic matrix around $x=0$. If $\mathrm{R}$ has only real eigenvalues, this package finds the intermediate matrices for the solutions using a power series expansion, around the singular point.
The general solutions, up to integration constants, are of the form
```math
 SolutionsMatrix(x) = TUMatrix.PMatrix.x^{SMatrix}\,.
```
The method can be found in chapter 4 of E. A. Coddington and N. Levinson, *Theory of ordinary differential equations* (McGraw-Hill, New York, 1955) or chapter 6.2 of E. A. Coddington an R. Carlson, *Linear Ordinary Differential Equations* (Society for Industrial and Applied Mathematics, Philadelphia, 1997).

A mathematical notebook with examples is provided. Some of the provided examples show how the method can be used to find the analytical eigenfunctions of the radial perturbations of static, spherically symmetric spacetimes in the theory of General Relativity.
