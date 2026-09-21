#pragma once
#include <iostream>
#include <iomanip>
#include "Classes2D/Mesh2D/Mesh2D.hpp"
#include "Classes2D/SolverSettings/SolutionParameters.hpp"
#include "Classes2D/Matrices/RHS.hpp"
#include <string>
#include <vector>
#include <cstddef>

enum class  Scheme
{
    Laasonen,
    CrankNicolson
};

class CoefficientMatrix
{
public:
    CoefficientMatrix(
    Scheme ActiveScheme,
    const Mesh2D& Mesh2D_Obj,
    const SolutionParameters<double> DiffNumberX_Obj,
    const SolutionParameters<double> DiffNumberY_Obj,
    double Default_Value_ = 0.0);

[[nodiscard]] std::size_t GetNumRows() const;   //[[nodiscard]] : The output value of this function should not be left unused.
[[nodiscard]] std::size_t GetNumColumns() const;
[[nodiscard]] double GetValue(std::size_t RowIndex, std::size_t ColumnIndex) const;

void SetValue(std::size_t RowIndex,std::size_t ColumnIndex,double Value_);

 void PrintFullSystem(const RHS& rhs) const;
 void SetBoundaryConditions();


private:
    std::size_t nrows;
    std::size_t ncols;
    std::size_t NN;
    std::vector<std::vector<double>> AValues;
    double rx;
    double ry;

};
