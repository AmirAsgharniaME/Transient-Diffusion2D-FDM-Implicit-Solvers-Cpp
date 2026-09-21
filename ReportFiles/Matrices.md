# Linear System Formulation and Dirichlet Boundary Conditions

In implicit finite difference formulations for the one-dimensional transient diffusion equation, the spatial coupling of unknown nodal values at the advanced time level $n+1$ leads to a coupled system of linear algebraic equations. This section presents the general coefficient matrix assembled for the interior nodes and establishes how Dirichlet boundary conditions are integrated into the algebraic framework.

## A. System of Equations for Interior Nodes

Consider a one-dimensional domain discretized into $N$ grid points indexed by $j = 0, 1, \dots, N-1$. The physical boundaries are located at $j = 0$ (lower wall) and $j = N-1$ (upper wall). Consequently, the interior computational domain consists of $N-2$ unknown nodes:

$$
j = 1, 2, \dots, N-2
$$

Evaluating an implicit discretization stencil at each interior node couples the unknown velocity $u_j^{n+1}$ to its adjacent spatial neighbors $u_{j-1}^{n+1}$ and $u_{j+1}^{n+1}$. For instance, employing the standard fully implicit (Laasonen) scheme yields the following discrete equation at node $j$:

$$
-d \, u_{j-1}^{n+1} + (1 + 2d) \, u_j^{n+1} - d \, u_{j+1}^{n+1} = u_j^n
$$ {#eq:implicit-interior-discrete}

where the dimensionless diffusion number $d$ is defined as:

$$
d = \frac{\nu \Delta t}{(\Delta y)^2}
$$ {#eq:implicit-diffusion-number}

Writing Equation ([-@eq:implicit-interior-discrete]) for all interior nodes $j = 1, 2, \dots, N-2$ simultaneously produces an $(N-2) \times (N-2)$ linear system of equations:

$$
\begin{bmatrix}
1+2d & -d    & 0     & \cdots & 0 \\
-d   & 1+2d  & -d    & \ddots & \vdots \\
0    & -d    & 1+2d  & \ddots & 0 \\
\vdots & \ddots & \ddots & \ddots & -d \\
0 & \cdots & 0 & -d & 1+2d
\end{bmatrix}
\begin{bmatrix}
u_1^{n+1} \\
u_2^{n+1} \\
u_3^{n+1} \\
\vdots \\
u_{N-2}^{n+1}
\end{bmatrix}
=
\begin{bmatrix}
u_1^n \\
u_2^n \\
u_3^n \\
\vdots \\
u_{N-2}^n
\end{bmatrix}
+
\begin{bmatrix}
d \, u_0^{n+1} \\
0 \\
\vdots \\
0 \\
d \, u_{N-1}^{n+1}
\end{bmatrix}
$$ {#eq:implicit-interior-system}

In Equation ([-@eq:implicit-interior-system]), the boundary contributions $u_0^{n+1}$ and $u_{N-1}^{n+1}$ appear directly on the right-hand side because their values are prescribed by the boundary conditions.


## B. Incorporation of Dirichlet Boundary Conditions

For the physical Couette flow problem, the boundary values at the lower and upper boundaries are strictly prescribed by constant Dirichlet boundary conditions:

$$
u_0^{n+1} = U_{\mathrm{bottom}}, \qquad u_{N-1}^{n+1} = U_{\mathrm{top}}
$$ {#eq:dirichlet-boundary-conditions}

There are two primary approaches for incorporating these boundary conditions into the coefficient matrix and algebraic system.

### Approach A: Reduced Interior System ($(N-2) \times (N-2)$)

In the reduced interior formulation, the known boundary values are excluded from the unknown solution vector. As shown in Equation ([-@eq:implicit-interior-system]), their contributions are absorbed directly into the right-hand side vector:

$$
\mathbf{A}_{\mathrm{int}} \, \mathbf{u}_{\mathrm{int}}^{n+1} = \mathbf{b}_{\mathrm{int}}^n
$$ {#eq:reduced-matrix-form}

where $\mathbf{u}_{\mathrm{int}}^{n+1} = [u_1^{n+1}, u_2^{n+1}, \dots, u_{N-2}^{n+1}]^T$ is of dimension $N-2$, and the right-hand side vector elements are given by:

$$
\begin{aligned}
b_1 &= u_1^n + d \, U_{\mathrm{bottom}} \\
b_j &= u_j^n \quad \text{for } j = 2, 3, \dots, N-3 \\
b_{N-2} &= u_{N-2}^n + d \, U_{\mathrm{top}}
\end{aligned}
$$ {#eq:reduced-rhs-terms}

This approach minimizes the algebraic system to $(N-2)$ equations, but it introduces an index offset between the algebraic vector index (from $0$ to $N-3$) and the physical spatial grid index (from $1$ to $N-2$).



### Approach B: Global Formulation ($N \times N$)

In the global system formulation, all computational nodes across the domain—including the boundary nodes—are preserved in the solution vector $\mathbf{u}^{n+1} \in \mathbb{R}^N$.

To close the system, the first and last rows of the coefficient matrix are formulated as identity equations corresponding to the Dirichlet boundary conditions in Equation ([-@eq:dirichlet-boundary-conditions]):

$$
\begin{aligned}
\text{Row } 0: \quad & 1 \cdot u_0^{n+1} + 0 \cdot u_1^{n+1} + \dots + 0 \cdot u_{N-1}^{n+1} = U_{\mathrm{bottom}} \\
\text{Row } N-1: \quad & 0 \cdot u_0^{n+1} + \dots + 0 \cdot u_{N-2}^{n+1} + 1 \cdot u_{N-1}^{n+1} = U_{\mathrm{top}}
\end{aligned}
$$ {#eq:global-boundary-rows}

Assembling the complete $N \times N$ global system yields:

$$
\begin{bmatrix}
1 & 0 & 0 & \cdots & 0 & 0 \\
-d & 1+2d & -d & \cdots & 0 & 0 \\
0 & -d & 1+2d & \ddots & \vdots & \vdots \\
\vdots & \ddots & \ddots & \ddots & -d & 0 \\
0 & \cdots & 0 & -d & 1+2d & -d \\
0 & 0 & \cdots & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
u_0^{n+1} \\
u_1^{n+1} \\
u_2^{n+1} \\
\vdots \\
u_{N-2}^{n+1} \\
u_{N-1}^{n+1}
\end{bmatrix}
=
\begin{bmatrix}
U_{\mathrm{bottom}} \\
u_1^n \\
u_2^n \\
\vdots \\
u_{N-2}^n \\
U_{\mathrm{top}}
\end{bmatrix}
$$ {#eq:global-linear-system}



### C. Selection for C++ Implementation

While Equation ([-@eq:reduced-matrix-form]) and Equation ([-@eq:global-linear-system]) are mathematically equivalent for all interior nodes, the global $N \times N$ formulation defined in Equation ([-@eq:global-linear-system]) is adopted in the current C++ solver architecture.

The primary advantages of this formulation are:

\begin{enumerate}
    \item \textbf{Uniform Indexing:} The row index of the algebraic system corresponds directly to the physical grid index $j$, which eliminates off-by-one indexing errors.
    \item \textbf{Simplified Memory Layout:} The solution vector contains the complete domain profile, eliminating the need to manually reinsert boundary values after each linear solve.
    \item \textbf{Boundary Condition Extensibility:} Alternative boundary conditions (such as Neumann or Robin boundaries) can be implemented simply by adjusting the discrete coefficients in row $0$ or row $N-1$, without altering the dimension or structural mapping of the linear solver.
\end{enumerate}

