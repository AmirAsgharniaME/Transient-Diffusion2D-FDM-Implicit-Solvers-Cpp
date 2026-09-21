#pragma once

#include <iostream>
#include <iomanip>
#include <cstddef>
#include <vector>

#include "Classes2D/Field2D/Field2D.hpp"
#include "Classes2D/SolverSettings/SolutionParameters.hpp"
#include "EssentialHeaders/SolverInputs.hpp"

class RHS
{


private:
    std::size_t nrows;
    std::size_t ncols;
    std::size_t NN;
    std::vector<double> RHSValues;

    double d_x; //x Difussion Number
    double d_y; //y Difussion Number
  


public:

    //RHS For Lassonen Method
    explicit RHS(
        const Field2D& Field2D_n_Obj,
        double Default_Initial_Value_ = 0.0
    );


    //RHS For CrankNicolson Method
    explicit RHS(
        const Field2D& Field2D_n_Obj,
        const SolutionParameters<double>& DiffNumerX_Obj,
        const SolutionParameters<double>& DiffNumerY_Obj,
        double Default_Initial_Value_ = 0.0      
    );

    
    /*****Public Setter*****/ 
    void SetValue(std::size_t PIndex, double Value_);
    void SetAllSingleValue(double Value_);

    /*****Public Getters*****/ 
     std::size_t GetSize() const;
     double GetValue(std::size_t PIndex) const;
     void SetBoundaryConditions();


    



};

void swap(Field2D& Left, Field2D& Right) noexcept;

