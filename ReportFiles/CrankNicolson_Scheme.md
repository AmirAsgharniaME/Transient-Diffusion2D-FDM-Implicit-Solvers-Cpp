# Crank-Nicolson Scheme

The Crank-Nicolson scheme is an implicit finite difference method that achieves second-order accuracy in time by evaluating the diffusion term at the arithmetic average of the current and future time levels. In contrast to the Laasonen scheme, which uses the diffusion term only at time level $n+1$, the Crank-Nicolson method blends the contributions from time levels $n$ and $n+1$.

For the one-dimensional transient diffusion equation, the Crank-Nicolson formulation can be written as the average of the spatial diffusion term at the two adjacent time levels:

$$
\left( \frac{\partial u}{\partial t} \right)_j^{n+\frac{1}{2}}
=
\frac{\nu}{2}
\left[
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n+1}
+
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n}
\right]
$$ {#eq:crank-nicolson-governing-form}

This midpoint form is the key distinction between Crank-Nicolson and Laasonen. Both methods are implicit, but Crank-Nicolson is more accurate in time.

## A. Time Derivative Discretization

To derive the finite difference approximation for the time derivative, expand $u_j^{n+1}$ and $u_j^n$ about the midpoint time level $t^{n+\frac{1}{2}}$.

Expanding $u_j^{n+1}$ about $t^{n+\frac{1}{2}}$ gives:

$$
u_j^{n+1}
=
u_j^{n+\frac{1}{2}}
+
\frac{\Delta t}{2}
\left(
\frac{\partial u}{\partial t}
\right)_j^{n+\frac{1}{2}}
+
\frac{(\Delta t)^2}{8}
\left(
\frac{\partial^2 u}{\partial t^2}
\right)_j^{n+\frac{1}{2}}
+
\mathcal{O}((\Delta t)^3)
$$ {#eq:crank-nicolson-time-taylor-forward}

Similarly, expanding $u_j^n$ about $t^{n+\frac{1}{2}}$ gives:

$$
u_j^{n}
=
u_j^{n+\frac{1}{2}}
-
\frac{\Delta t}{2}
\left(
\frac{\partial u}{\partial t}
\right)_j^{n+\frac{1}{2}}
+
\frac{(\Delta t)^2}{8}
\left(
\frac{\partial^2 u}{\partial t^2}
\right)_j^{n+\frac{1}{2}}
-
\mathcal{O}((\Delta t)^3)
$$ {#eq:crank-nicolson-time-taylor-backward}

Subtracting Equation ([-@eq:crank-nicolson-time-taylor-backward]) from Equation ([-@eq:crank-nicolson-time-taylor-forward]) yields:

$$
u_j^{n+1} - u_j^n
=
\Delta t
\left(
\frac{\partial u}{\partial t}
\right)_j^{n+\frac{1}{2}}
+
\mathcal{O}((\Delta t)^3)
$$ {#eq:crank-nicolson-time-subtraction}

Therefore, the midpoint time derivative can be approximated as:

$$
\left(
\frac{\partial u}{\partial t}
\right)_j^{n+\frac{1}{2}}
=
\frac{u_j^{n+1} - u_j^n}{\Delta t}
+
\mathcal{O}((\Delta t)^2)
$$ {#eq:crank-nicolson-time-derivative}

Thus, the Crank-Nicolson scheme is second-order accurate in time, whereas Laasonen is only first-order accurate in time.

## B. Spatial Derivative Discretization

For the diffusion term, the second derivative is approximated at both time levels $n$ and $n+1$ using the central difference formula.

At time level $n+1$:

$$
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n+1}
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
\mathcal{O}((\Delta y)^2)
$$ {#eq:crank-nicolson-spatial-nplus1}

At time level $n$:

$$
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n}
=
\frac{
u_{j+1}^{n}
-
2u_j^{n}
+
u_{j-1}^{n}
}{
(\Delta y)^2
}
+
\mathcal{O}((\Delta y)^2)
$$ {#eq:crank-nicolson-spatial-n}

Taking the arithmetic average of Equations ([-@eq:crank-nicolson-spatial-nplus1]) and ([-@eq:crank-nicolson-spatial-n]) gives:

$$
\frac{1}{2}
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n+1}
+
\frac{1}{2}
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n}
=
\frac{1}{2(\Delta y)^2}
\left[
u_{j+1}^{n+1}
-
2u_j^{n+1}
+
u_{j-1}^{n+1}
+
u_{j+1}^{n}
-
2u_j^{n}
+
u_{j-1}^{n}
\right]
+
\mathcal{O}((\Delta y)^2)
$$ {#eq:crank-nicolson-spatial-average}

This confirms that the Crank-Nicolson scheme retains second-order spatial accuracy, the same as Laasonen and FTCS.

## C. Discrete Algebraic Equation

Substituting Equation ([-@eq:crank-nicolson-time-derivative]) and Equation ([-@eq:crank-nicolson-spatial-average]) into Equation ([-@eq:crank-nicolson-governing-form]) gives:

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
\frac{1}{2(\Delta y)^2}
\left[
u_{j+1}^{n+1}
-
2u_j^{n+1}
+
u_{j-1}^{n+1}
+
u_{j+1}^{n}
-
2u_j^{n}
+
u_{j-1}^{n}
\right]
$$ {#eq:crank-nicolson-discrete-form}

Defining the diffusion number as:

$$
d
=
\frac{\nu \Delta t}{(\Delta y)^2}
$$ {#eq:crank-nicolson-diffusion-number}

the Crank-Nicolson scheme can be rearranged into the standard implicit form:

$$
-\frac{d}{2}u_{j-1}^{n+1}
+
(1+d)u_j^{n+1}
-
\frac{d}{2}u_{j+1}^{n+1}
=
\frac{d}{2}u_{j-1}^{n}
+
(1-d)u_j^{n}
+
\frac{d}{2}u_{j+1}^{n}
$$ {#eq:crank-nicolson-final-form}

Equation ([-@eq:crank-nicolson-final-form]) is the discrete algebraic equation for Crank-Nicolson. It has the same tridiagonal coupling pattern as Laasonen, but the right-hand side now also contains information from time level $n$.



## D. Accuracy of the Crank-Nicolson Scheme

The temporal truncation error is:

$$
\mathcal{O}((\Delta t)^2)
$$

The spatial truncation error is:

$$
\mathcal{O}((\Delta y)^2)
$$

Therefore, the Crank-Nicolson scheme is:

- second-order accurate in time;
- second-order accurate in space.

The combined local truncation error is:

$$
\mathcal{T}_{\mathrm{CN}}
=
\mathcal{O}((\Delta t)^2)
+
\mathcal{O}((\Delta y)^2)
$$ {#eq:crank-nicolson-local-error}

This is the main accuracy advantage of Crank-Nicolson over Laasonen, which is only first-order accurate in time.


## E. Comparison Between Crank-Nicolson and Laasonen

The two implicit schemes used in this report are Laasonen and Crank-Nicolson. Both methods lead to tridiagonal linear systems and both are unconditionally stable for the one-dimensional diffusion equation. Their principal difference lies in how the diffusion term is treated in time.

### Laasonen
- first-order accurate in time;
- second-order accurate in space;
- fully implicit at time level $n+1$;
- more strongly damping;
- numerically robust and monotone for diffusion-dominated problems.

### Crank-Nicolson
- second-order accurate in time;
- second-order accurate in space;
- averages the diffusion term between $n$ and $n+1$;
- less numerical damping than Laasonen;
- may show oscillations for large $d$.

In summary, Crank-Nicolson is more accurate than Laasonen, but Laasonen is often more diffusive and therefore more forgiving for strongly dissipative transient problems. If the goal is higher temporal accuracy, Crank-Nicolson is preferable. If the goal is stronger numerical smoothing and simpler behavior, Laasonen can be advantageous.

## I. Comparison of Implicit and Explicit Schemes

The numerical schemes discussed in this report can be divided into explicit and implicit classes.

### Explicit schemes
- FTCS
- Richardson
- DuFort-Frankel

### Implicit schemes
- Laasonen
- Crank-Nicolson

The main differences are:

| Property | Explicit schemes | Implicit schemes |
|---|---|---|
| Unknowns at new time level | Computed directly | Found by solving a linear system |
| Stability | Usually restricted or problematic | Typically unconditionally stable for diffusion |
| Computational cost per step | Low | Higher |
| Time-step flexibility | Limited by stability | Much larger |
| Implementation | Simpler | More involved |
| Accuracy control | Often lower or unstable | Can be higher and more robust |

For the present diffusion problem:

- FTCS is simple but conditionally stable.
- Richardson is formally attractive but unusable because it is unconditionally unstable.
- DuFort-Frankel is explicit and unconditionally stable, but it uses three time levels and requires a startup procedure.
- Laasonen is robust and unconditionally stable, but only first-order accurate in time.
- Crank-Nicolson is the most accurate among the schemes considered here, with second-order accuracy in both time and space, but it may require care when large time steps create oscillatory behavior.

Therefore, the implicit schemes provide better stability properties, while the explicit schemes provide simpler updates but weaker stability guarantees.

## J. Summary of the Crank-Nicolson Scheme

The main properties of the Crank-Nicolson scheme are:

- It is an implicit finite difference method.
- It evaluates the diffusion term at the average of time levels $n$ and $n+1$.
- It is second-order accurate in time.
- It is second-order accurate in space.
- It produces a tridiagonal linear system at every time step.
- It is unconditionally stable for the one-dimensional diffusion equation.
- It may produce oscillations for large diffusion numbers.
- It is more accurate in time than Laasonen.
- It requires a tridiagonal solver such as the Thomas algorithm.

Therefore, Crank-Nicolson provides a high-accuracy implicit alternative for the transient Plane Couette Flow problem, especially when accuracy is more important than the extra damping provided by Laasonen.
