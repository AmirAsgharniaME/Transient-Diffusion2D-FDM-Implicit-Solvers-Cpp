#include "Classes2D/Matrices/CoefficientMatrix.hpp"
#include <cstddef>

CoefficientMatrix::CoefficientMatrix(
    Scheme ActiveScheme,
    const Mesh2D& Mesh2D_Obj,
    const SolutionParameters<double> DiffNumberX_Obj,
    const SolutionParameters<double> DiffNumberY_Obj,
    double Default_Value_)
    :nrows(Mesh2D_Obj.GetNumOfYNodes()),
     ncols(Mesh2D_Obj.GetNumOfXNodes()),
     NN(nrows*ncols),
     AValues(NN, std::vector<double>(NN, Default_Value_)),
     rx(DiffNumberX_Obj.GetValue()),
     ry(DiffNumberY_Obj.GetValue())
    {

    if (ActiveScheme == Scheme::Laasonen)
    {
    // ================
    // Laasonen Method
    // ================
    // EQ : (1 + 2*dx + 2*dy) * u_(n+1)[j][i]
    //     - dx * u_(n+1)[j][i - 1]
    //     - dx * u_(n+1)[j][i + 1]
    //     - dy * u_(n+1)[j - 1][i]
    //     - dy * u_(n+1)[j + 1][i]
    //     = u_(n)[j][i]
    // 1D numbering : p = j * Nx + i
    // RHS = U[j][i][n]
    // Create A for Laasonen
    // ========================================
    std::size_t Nx=ncols;
    std::size_t Ny=nrows;
    const double dx = rx;
    const double dy = ry;
    // Iterate through the interior nodes
    for (std::size_t j = 1; j < Ny - 1; ++j)
    {
        for (std::size_t i = 1; i < Nx - 1; ++i)
        {
            const std::size_t p = j * Nx + i;
            AValues[p][p - 1] =  -dx;                         // West Neighbor (i - 1)
            AValues[p][p] =      1.0 + 2.0 * dx + 2.0 * dy;   // Center (i, j)
            AValues[p][p + 1] =  -dx;                         // East Neighbor (i + 1)
            AValues[p][p - Nx] = -dy;                         // South / Bottom Neighbor (j - 1)
            AValues[p][p + Nx] = -dy;                         // North / Top Neighbor (j + 1)
    
        }
        
    } 


        }


    // ===============================================
    // 2D Diffusion Equation - Crank-Nicolson Method
    // ==============================================
//  * Full equation:
//  *   -0.5*dx*u[j][i-1][n+1] - 0.5*dx*u[j][i+1][n+1] + (1+dx+dy)*u[j][i][n+1]
//  *   -0.5*dy*u[j-1][i][n+1] - 0.5*dy*u[j+1][i][n+1]
//  *   =
//  *   0.5*dx*u[j][i-1][n] + 0.5*dx*u[j][i+1][n] + (1-dx-dy)*u[j][i][n]
//  *   +0.5*dy*u[j-1][i][n] + 0.5*dy*u[j+1][i][n]


    if (ActiveScheme == Scheme::CrankNicolson)
    {
    std::size_t Nx=ncols;
    std::size_t Ny=nrows;
    const double dx = rx;
    const double dy = ry;
    // Iterate through the interior nodes
    for (std::size_t j = 1; j < Ny - 1; ++j)
    {
        for (std::size_t i = 1; i < Nx - 1; ++i)
        {
            const std::size_t p = j * Nx + i;
            AValues[p][p - 1] =  -0.5*dx;                 // West Neighbor (i - 1)
            AValues[p][p] =      1.0 + dx + dy;          // Center (i, j)
            AValues[p][p + 1] =  - 0.5*dx;               // East Neighbor (i + 1)
            AValues[p][p - Nx] = -0.5*dy;                 // South / Bottom Neighbor (j - 1)
            AValues[p][p + Nx] = -0.5*dy;                  // North / Top Neighbor (j + 1)
    
        }

        }

    }

    }

std::size_t CoefficientMatrix::GetNumRows() const
{
    // size of a 2D vector is num of rows.
    return AValues.size();
}

std::size_t CoefficientMatrix::GetNumColumns() const
{
    // size of first row (1D vector) in num of columns.
    // Equal to AValues[0].size();
    return AValues.empty() ? 0 : AValues.front().size(); 
}

double CoefficientMatrix::GetValue(std::size_t RowIndex,std::size_t ColumnIndex) const
{
    return AValues[RowIndex][ColumnIndex];
}

void CoefficientMatrix::SetValue(
    std::size_t RowIndex,
    std::size_t ColumnIndex,
    double Value_)
{
    AValues[RowIndex][ColumnIndex] = Value_;
}

void CoefficientMatrix::PrintFullSystem(const RHS& rhs) const
{
    // Print column header with node coordinates, location, RHS value, and matrix row entries
    std::cout << "Node(j,i) | Location |   RHS Value   | Matrix Row Content" << std::endl;
    std::cout << "----------------------------------------------------------------------------------------" << std::endl;

    for (std::size_t row = 0; row < NN; row++)
    {
        // Convert the flat matrix row index back to 2D grid coordinates (j, i)
        std::size_t j = row / ncols; // nrows is Ny, ncols is Nx
        std::size_t i = row % ncols;

        // Determine boundary or interior location based on defined priority: Top, Bottom, Right, Left
        std::string location = "Interior";

        if (j == nrows - 1)
        {
            location = "Top     ";
        }
        else if (j == 0)
        {
            location = "Bottom  ";
        }
        else if (i == ncols - 1)
        {
            location = "Right   ";
        }
        else if (i == 0)
        {
            location = "Left    ";
        }

        // Print node coordinates, location label, and the corresponding RHS value
        std::cout << "(" << j << "," << i << ") | " 
                  << location << " | "
                  << std::scientific << std::setprecision(2) << std::setw(13) << rhs.GetValue(row) << " | ";

        // Print the row elements of the coefficient matrix
        for (std::size_t col = 0; col < NN; col++)
        {
            if (std::abs(AValues[row][col]) < 1e-12)
            {
                std::cout << "       .      ";
            }
            else
            {
                std::cout << std::scientific
                          << std::setprecision(2)
                          << std::setw(12)
                          << AValues[row][col] << " ";
            }
        }
        std::cout << std::endl;
    }
}





void CoefficientMatrix::SetBoundaryConditions()
{

    std::size_t Nx=ncols;
    std::size_t Ny=nrows;

    // 1. Top Boundary (j = Ny - 1)
    {
        const std::size_t j = Ny - 1;
        for (std::size_t i = 0; i < Nx; ++i)
        {
            const std::size_t p = j * Nx + i;
            AValues[p][p] = 1.0;
        }
    }

    // 2. Bottom Boundary (j = 0)
    {
        const std::size_t j = 0;
        for (std::size_t i = 0; i < Nx; ++i)
        {
            const std::size_t p = j * Nx + i;
            AValues[p][p] = 1.0;
        }
    }

    // 3. Right Boundary (i = Nx - 1)
    {
        const std::size_t i = Nx - 1;
        for (std::size_t j = 0; j < Ny; ++j)
        {
            const std::size_t p = j * Nx + i;
            AValues[p][p] = 1.0;
        }
    }

    // 4. Left Boundary (i = 0)
    {
        const std::size_t i = 0;
        for (std::size_t j = 0; j < Ny; ++j)
        {
            const std::size_t p = j * Nx + i;
            AValues[p][p] = 1.0;
        }
    }
}