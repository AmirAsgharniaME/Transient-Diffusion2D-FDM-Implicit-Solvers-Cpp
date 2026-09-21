#include "Classes2D/Matrices/RHS.hpp"
//#include <algorithm>


    //RHS For Lassonen Method
    RHS::RHS(
        const Field2D& Field2D_n_Obj,
        double Default_Initial_Value_ )
        :nrows(Field2D_n_Obj.GetNrows()),
         ncols(Field2D_n_Obj.GetNcols()),
         NN(nrows*ncols),
         RHSValues(NN, Default_Initial_Value_)
    {
        // ==============
        // RHS For Lassonen Method
        // ===============
        // EQ :-d_x * u(j-1, i)^{n+1} - d_y * u(j, i-1)^{n+1} + (1 + 2*d_x + 2*d_y) * u(j, i)^{n+1} - d_x * u(j+1, i)^{n+1} - d_y * u(j, i+1)^{n+1} = u(j, i)^n
        // RHS = U[j][i][n]
        // Create RHS for Laasonen Method
           std::size_t p;
           std::size_t Nx = ncols;

            // for internal nodes
            for (std::size_t j = 1; j < nrows-1; j++)
            {
                for (std::size_t i = 1; i < ncols-1; i++)
                {
                    p = j*Nx + i;
                    RHSValues[p] = Field2D_n_Obj.GetValue(j,i);
                    
                }
                
            }

            //Top Wall
            //for (std::size_t i = 1; i < ncols-1; i++)

    }


    //RHS For CrankNicolson Method
     RHS::RHS(
        const Field2D& Field2D_n_Obj,
        const SolutionParameters<double>& DiffNumerX_Obj,
        const SolutionParameters<double>& DiffNumerY_Obj,
        double Default_Initial_Value_     
    )
    :nrows(Field2D_n_Obj.GetNrows()),
     ncols(Field2D_n_Obj.GetNcols()),
     NN(nrows*ncols),
     RHSValues(NN, Default_Initial_Value_),
     d_x(DiffNumerX_Obj.GetValue()),
     d_y(DiffNumerY_Obj.GetValue())

     {
        // ===================
        // CrankNicolson Method
        // ====================

        // LHS (Unknowns at time step n+1):
        // -0.5*d_x * u(j, i-1)^{n+1} - 0.5*d_y * u(j-1, i)^{n+1} + (1 + d_x + d_y) * u(j, i)^{n+1} - 0.5*d_x * u(j, i+1)^{n+1} - 0.5*d_y * u(j+1, i)^{n+1}

        // RHS (Knowns at time step n):
        // 0.5*d_x * u(j, i-1)^n + 0.5*d_y * u(j-1, i)^n + (1 - d_x - d_y) * u(j, i)^n + 0.5*d_x * u(j, i+1)^n + 0.5*d_y * u(j+1, i)^n
        //=============================================================================================================================

        // Create RHS for CrankNicolson Method
           std::size_t p;
           std::size_t Nx = ncols;

            // for internal nodes
            for (std::size_t j = 1; j < nrows-1; j++)
            {
                for (std::size_t i = 1; i < ncols-1; i++)
                { 
                    p = j*Nx + i;

                   // for internal nodes
                   double value = 0.5*d_x * Field2D_n_Obj.GetValue(j, i-1)
                   + 0.5*d_y * Field2D_n_Obj.GetValue(j-1, i) 
                   + (1 - d_x - d_y) * Field2D_n_Obj.GetValue(j, i) 
                   + 0.5*d_x * Field2D_n_Obj.GetValue(j, i+1) 
                   + 0.5*d_y * Field2D_n_Obj.GetValue(j+1, i);
                  
                   RHSValues[p] = value;
                   
                }
                
            }

     }


 void RHS::SetValue(std::size_t PIndex, double Value_)
    {
        RHSValues[PIndex] = Value_;
    }
    
double RHS::GetValue(std::size_t PIndex) const
{

    return RHSValues[PIndex];
}

std::size_t RHS::GetSize() const
{
return RHSValues.size();
}

void  RHS::SetBoundaryConditions()
{
    std::size_t Nx=ncols;
    std::size_t Ny=nrows;
    std::size_t p;

    // 1. Top Boundary (j = Ny - 1)
    {
        const std::size_t j = Ny - 1;
        for (std::size_t i = 0; i < Nx; ++i)
        {
            p = j*Nx + i;
            RHSValues[p] = SolverInputs::TopWall_Temp_Values[i];
        }
    }
    // 2. Bottom Boundary (j = 0)
    {
        const std::size_t j = 0;
        for (std::size_t i = 0; i < Nx; ++i)
        {
            p = j*Nx + i;
            RHSValues[p] = SolverInputs::BottomWall_Temp_Values[i];
        }
    }

    // 3. Right Boundary (i = Nx - 1)
    {
        const std::size_t i = Nx - 1;
        for (std::size_t j = 0; j < Ny; ++j)
        {
            p = j*Nx + i;
            RHSValues[p] = SolverInputs::RightWall_temp_Values[j];
        }
    }

    // 4. Left Boundary (i = 0)
    {
        const std::size_t i = 0;
        for (std::size_t j = 0; j < Ny; ++j)
        {
            p = j*Nx + i;
            RHSValues[p] = SolverInputs::LeftWall_Temp_Values[j];
        }
    }

    
}