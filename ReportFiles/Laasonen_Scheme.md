# Laasonen Scheme

The Implicit Euler scheme, commonly known as the Laasonen scheme, is an implicit finite difference method for solving the one-dimensional transient diffusion equation.

Unlike the FTCS and Richardson schemes, which evaluate the spatial derivative using only previously known values, the Laasonen scheme evaluates the spatial derivative at the new time level $n+1$. Consequently, the unknown values at the new time level appear simultaneously in the discretized equation, and a system of linear algebraic equations must be solved at every time step.

The governing one-dimensional diffusion equation, evaluated at the grid point $(y_j,t^{n+1})$, gives:

$$
\left(\frac{\partial u}{\partial t}\right)_j^{n+1}
=
\nu
\left(\frac{\partial^2 u}{\partial y^2}\right)_j^{n+1}
$$ {#eq:laasonen-discrete-notation}

The Laasonen scheme uses:

- a first-order backward difference in time;
- a second-order central difference in space;
- an implicit evaluation of the diffusion term at time level $n+1$;
- a tridiagonal linear system at each time step.

## A. Discretization of the First-Order Temporal Derivative

To approximate the first-order temporal derivative at time level $n+1$, the value $u_j^n$ is expanded about the point $(y_j,t^{n+1})$:

$$
u_j^n
=
u_j^{n+1}
-
\left(\frac{\partial u}{\partial t}\right)_j^{n+1}\Delta t
+
\frac{1}{2!}
\left(\frac{\partial^2 u}{\partial t^2}\right)_j^{n+1}
(\Delta t)^2
-
\frac{1}{3!}
\left(\frac{\partial^3 u}{\partial t^3}\right)_j^{n+1}
(\Delta t)^3
+
\mathcal{O}\left((\Delta t)^4\right)
$$ {#eq:laasonen-backward-temporal-taylor}

Rearranging Equation ([-@eq:laasonen-backward-temporal-taylor]) to isolate the first-order temporal derivative gives:

$$
\left(\frac{\partial u}{\partial t}\right)_j^{n+1}
=
\frac{
u_j^{n+1}
-
u_j^n
}{
\Delta t
}
+
\mathcal{O}(\Delta t)
$$ {#eq:laasonen-backward-time}

Therefore, the backward difference approximation is first-order accurate in time.

## B. Discretization of the Second-Order Spatial Derivative

The second-order spatial derivative is evaluated at time level $n+1$. Expanding the neighboring values $u_{j+1}^{n+1}$ and $u_{j-1}^{n+1}$ about $(y_j,t^{n+1})$ gives:

$$
u_{j+1}^{n+1}
=
u_j^{n+1}
+
\left(\frac{\partial u}{\partial y}\right)_j^{n+1}\Delta y
+
\frac{1}{2!}
\left(\frac{\partial^2 u}{\partial y^2}\right)_j^{n+1}
(\Delta y)^2
+
\frac{1}{3!}
\left(\frac{\partial^3 u}{\partial y^3}\right)_j^{n+1}
(\Delta y)^3
+
\frac{1}{4!}
\left(\frac{\partial^4 u}{\partial y^4}\right)_j^{n+1}
(\Delta y)^4
+
\cdots
$$ {#eq:laasonen-forward-spatial-taylor}

$$
u_{j-1}^{n+1}
=
u_j^{n+1}
-
\left(\frac{\partial u}{\partial y}\right)_j^{n+1}\Delta y
+
\frac{1}{2!}
\left(\frac{\partial^2 u}{\partial y^2}\right)_j^{n+1}
(\Delta y)^2
-
\frac{1}{3!}
\left(\frac{\partial^3 u}{\partial y^3}\right)_j^{n+1}
(\Delta y)^3
+
\frac{1}{4!}
\left(\frac{\partial^4 u}{\partial y^4}\right)_j^{n+1}
(\Delta y)^4
+
\cdots
$$ {#eq:laasonen-backward-spatial-taylor}

Adding Equations ([-@eq:laasonen-forward-spatial-taylor]) and ([-@eq:laasonen-backward-spatial-taylor]) eliminates the odd-order spatial derivative terms:

$$
u_{j+1}^{n+1}
+
u_{j-1}^{n+1}
=
2u_j^{n+1}
+
\left(\frac{\partial^2 u}{\partial y^2}\right)_j^{n+1}
(\Delta y)^2
+
\frac{1}{12}
\left(\frac{\partial^4 u}{\partial y^4}\right)_j^{n+1}
(\Delta y)^4
+
\cdots
$$ {#eq:laasonen-central-taylor-sum}

Rearranging Equation ([-@eq:laasonen-central-taylor-sum]) gives the second-order central difference approximation:

$$
\left(\frac{\partial^2 u}{\partial y^2}\right)_j^{n+1}
=
\frac{
u_{j+1}^{n+1}
-
2u_j^{n+1}
+
u_{j-1}^{n+1}
}{
(\Delta y)^2
}
+
\mathcal{O}\left((\Delta y)^2\right)
$$ {#eq:laasonen-central-space}

Thus, the spatial discretization is second-order accurate.

## C. Discrete Algebraic Equation

Substituting the temporal approximation from Equation ([-@eq:laasonen-backward-time]) and the spatial approximation from Equation ([-@eq:laasonen-central-space]) into the governing equation in Equation ([-@eq:laasonen-discrete-notation]) gives:

$$
\frac{
u_j^{n+1}
-
u_j^n
}{
\Delta t
}
=
\nu
\left(
\frac{
u_{j+1}^{n+1}
-
2u_j^{n+1}
+
u_{j-1}^{n+1}
}{
(\Delta y)^2
}
\right)
$$ {#eq:laasonen-discrete-substitution}

Define the dimensionless diffusion number as:

$$
d=\frac{\nu\Delta t}{(\Delta y)^2}
$$ {#eq:laasonen-diffusion-number}

Multiplying Equation ([-@eq:laasonen-discrete-substitution]) by $\Delta t$ and using Equation ([-@eq:laasonen-diffusion-number]) yields:

$$
u_j^{n+1}
-
u_j^n
=
d
\left(
u_{j+1}^{n+1}
-
2u_j^{n+1}
+
u_{j-1}^{n+1}
\right)
$$ {#eq:laasonen-dimensionless-form}

Rearranging Equation ([-@eq:laasonen-dimensionless-form]) so that all unknown values at time level $n+1$ are placed on the left-hand side gives:

$$
-d\,u_{j-1}^{n+1}
+
(1+2d)\,u_j^{n+1}
-
d\,u_{j+1}^{n+1}
=
u_j^n
$$ {#eq:laasonen-tridiagonal-relation}

Equation ([-@eq:laasonen-tridiagonal-relation]) is the discrete algebraic form of the Laasonen scheme. It must be solved simultaneously for all spatial nodes at every time step. At every time step, the complete program must:

- Initialize the solution vector $u_j^n$ from the initial condition.
- Apply the boundary conditions on both sides.
- Assemble and solve the tridiagonal system.
- Update the solution to the next time level.

## D. Accuracy of the Laasonen Scheme

The temporal truncation error is obtained from Equation ([-@eq:laasonen-backward-time]):

$$
\text{Temporal truncation error}
=
\mathcal{O}(\Delta t)
$$

The spatial truncation error is obtained from Equation ([-@eq:laasonen-central-space]):

$$
\text{Spatial truncation error}
=
\mathcal{O}\left((\Delta y)^2\right)
$$

Therefore, the Laasonen scheme is:

- first-order accurate in time;
- second-order accurate in space.

The combined local truncation error can be written as:

$$
\mathcal{T}_{\mathrm{Laasonen}}
=
\mathcal{O}(\Delta t)
+
\mathcal{O}\left((\Delta y)^2\right)
$$ {#eq:laasonen-local-truncation-error}

The spatial accuracy of Laasonen is the same as that of FTCS, while its temporal accuracy is the same order as FTCS. The main difference between the two methods is their stability behavior and algebraic structure.





## E. Comparison with the Explicit Schemes

The main characteristics of the four schemes considered so far are summarized below:

| Scheme | Time discretization | Spatial discretization | Algebraic system | Temporal accuracy | Stability |
|---|---|---|---|---|---|
| FTCS | Forward explicit | Central | Not required | First-order | Conditionally stable, $d\leq\frac{1}{2}$ |
| Richardson | Central explicit | Central | Not required | Second-order | Unconditionally unstable |
| DuFort-Frankel | Central explicit | Modified central | Not required | Nominally second-order | Unconditionally stable |
| Laasonen | Backward implicit | Central | Tridiagonal system | First-order | Unconditionally stable |

The Laasonen scheme has the same formal temporal accuracy as the FTCS scheme, but it does not impose the FTCS stability restriction. Its main computational requirement is the solution of a tridiagonal linear system at every time step.

In contrast:

- FTCS is simple to implement because it requires only a direct explicit update.
- Richardson has second-order formal accuracy in both time and space, but it is unusable because of the continuous instability.
- DuFort-Frankel is explicit and unconditionally stable, but it uses three time levels and requires a starting procedure.
- Laasonen is implicit and unconditionally stable, but its solution requires a tridiagonal-system solver.

For the present transient Plane Couette Flow problem, Laasonen is a suitable method when stability and robustness are more important than minimizing the computational work per time step. Since its coefficient matrix is tridiagonal, the resulting system can be solved efficiently using a specialized tridiagonal solver, which will be described separately in a separately in a later section.

## F. Summary of the Laasonen scheme :

- It uses a first-order backward difference in time.
- It uses a second-order central difference in space.
- It is implicit because the unknown values at the new time level appear in the diffusion term.
- It produces a tridiagonal system at each step.
- It is first-order accurate in time and second-order accurate in space.
- It is unconditionally stable for the one-dimensional diffusion equation.
- It does not require a special initialization beyond the prescribed initial condition.
- It permits larger time steps than FTCS without violating the FTCS stability condition.
- Its linear system can be solved efficiently using a specialized tridiagonal solver.

Therefore, the Laasonen scheme provides a robust implicit alternative to the explicit FTCS, Richardson, and DuFort-Frankel schemes.
