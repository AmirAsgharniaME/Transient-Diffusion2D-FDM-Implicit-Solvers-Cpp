# Numerical Schemes

In this section, the transition from the continuous partial differential equation to discrete algebraic forms suitable for C++ implementation is presented.

To solve the one-dimensional transient diffusion equation, previously introduced in Theoretical Background as Equation ([-@eq:diffusion-final]):

$$
\frac{\partial u}{\partial t}
=
\nu \frac{\partial^2 u}{\partial y^2}
$$ {#eq:diffusion-final}

the continuous temporal and spatial derivatives must be approximated using finite difference expressions. Implicit numerical schemes are therefore introduced to obtain discrete formulations of the governing equation.

## Implicit Numerical Schemes

1. **Laasonen Scheme**
2. **CrankNicolson Scheme**
