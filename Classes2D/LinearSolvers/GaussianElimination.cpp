#include "Classes2D/LinearSolvers/GaussianElimination.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

void GaussianElimination::Solve_nPlus1(
    CoefficientMatrix& A,
    RHS& RHS,
    Field2D& Solution_nPlus1)
{
    const std::size_t NumARows = A.GetNumRows();
    const std::size_t NumACols = A.GetNumColumns();

    if (NumARows != NumACols)
    {
        throw std::invalid_argument("Solve: matrix A must be square.");
    }

    if (RHS.GetSize() != NumARows || Solution_nPlus1.GetNrows() * Solution_nPlus1.GetNcols() != NumARows)
    {
        throw std::invalid_argument("Solve: size of RHS and Solution must match matrix size.");
    }

    ForwardElimination(A, RHS);
    BackSubstitution(A, RHS, Solution_nPlus1);
}

void GaussianElimination::ForwardElimination(
    CoefficientMatrix& A,
    RHS& RHS)
{
    const std::size_t NumARows = A.GetNumRows();
    const std::size_t NumACols = A.GetNumColumns();

    if (NumARows != NumACols)
    {
        throw std::invalid_argument("ForwardElimination: matrix A must be square.");
    }

    if (RHS.GetSize() != NumARows)
    {
        throw std::invalid_argument("ForwardElimination: size of RHS must match matrix size.");
    }

    if (NumARows == 0)
    {
        return;
    }

    for (std::size_t PivotColIndex = 0; PivotColIndex < NumACols - 1; ++PivotColIndex)
    {
        std::size_t PivotRowIndex = PivotColIndex;

        const std::size_t Candidate_PivotRowIndex = FindPivotRowIndex(A, PivotColIndex);

        if (Candidate_PivotRowIndex != PivotRowIndex)
        {
            SwapRows(A, RHS, PivotRowIndex, Candidate_PivotRowIndex);
        }

        const double PivotValue = A.GetValue(PivotRowIndex, PivotColIndex);

        if (std::fabs(PivotValue) < 1.0e-14)
        {
            throw std::runtime_error("ForwardElimination: singular or nearly singular matrix detected.");
        }


        for (std::size_t RowIndex = PivotRowIndex + 1; RowIndex < NumARows; ++RowIndex)
        {
            const double Factor = A.GetValue(RowIndex, PivotColIndex) / PivotValue;

            A.SetValue(RowIndex, PivotColIndex, 0.0);

            for (std::size_t ColIndex = PivotColIndex + 1; ColIndex < NumACols; ++ColIndex)
            {
                const double NewValue =
                    A.GetValue(RowIndex, ColIndex) - Factor * A.GetValue(PivotRowIndex, ColIndex);
                A.SetValue(RowIndex, ColIndex, NewValue);
            }

            const double NewRHS =
                RHS.GetValue(RowIndex) - Factor * RHS.GetValue(PivotRowIndex);
            RHS.SetValue(RowIndex, NewRHS);
        }
    }

    if (std::fabs(A.GetValue(NumARows - 1, NumACols - 1)) < 1.0e-14)
    {
        throw std::runtime_error("ForwardElimination: zero pivot found on last row.");
    }
}

// Mathematical formulation of back-substitution:
// x_i = ( b_i - \sum_{j=i+1}^{N-1} A_{i,j} * x_j ) / A_{i,i}
//
// Where:
//   i         : Current row and equation index (from N-1 down to 0)
//   j         : Column indices to the right of the main diagonal (j > i)
//   A_{i,i}   : Pivot/diagonal entry of the current row
//   b_i       : Corresponding element in the right-hand side (RHS) vector
//   sum_upper : Cumulative sum of known upper terms: \sum_{j=i+1}^{N-1} A_{i,j} * x_j
//   x_i       : Computed unknown value stored in the solution vector

void GaussianElimination::BackSubstitution(
    const CoefficientMatrix& A,
    const RHS& RHS,
    Field2D& Solution_nPlus1)
{
    const std::size_t NumARows = A.GetNumRows();
    const std::size_t NumACols = A.GetNumColumns();

    if (NumARows != NumACols)
    {
        throw std::invalid_argument("BackSubstitution: matrix A must be square.");
    }

    if (RHS.GetSize() != NumARows || Solution_nPlus1.GetNrows() * Solution_nPlus1.GetNcols() != NumARows)
    {
        throw std::invalid_argument("BackSubstitution: size of RHS and Solution must match matrix size.");
    }

    const std::size_t ncols = Solution_nPlus1.GetNcols();

    // Mathematical formula:
    // x_i = ( b_i - \sum_{j=i+1}^{N-1} A_{i,j} * x_j ) / A_{i,i}

    // Iterate backward from the last row (N-1) down to the first row (0)
    for (std::ptrdiff_t row = static_cast<std::ptrdiff_t>(NumARows) - 1; row >= 0; --row)
    {
        const std::size_t i = static_cast<std::size_t>(row);

        // Calculate the sum of known upper terms: \sum_{j=i+1}^{N-1} A_{i,j} * x_j
        double sum_upper = 0.0;
        for (std::size_t j = i + 1; j < NumACols; ++j)
        {
            // Map 1D column index j to 2D grid coordinates (sol_j, sol_i)
            const std::size_t sol_j = j / ncols;
            const std::size_t sol_i = j % ncols;
            sum_upper += A.GetValue(i, j) * Solution_nPlus1.GetValue(sol_j, sol_i);
        }

        const double A_ii = A.GetValue(i, i);

        if (std::fabs(A_ii) < 1.0e-14)
        {
            throw std::runtime_error("BackSubstitution: zero diagonal entry detected.");
        }

        const double b_i = RHS.GetValue(i);

        // Map 1D row index i to 2D grid coordinates (target_j, target_i)
        const std::size_t target_j = i / ncols;
        const std::size_t target_i = i % ncols;

        // x_i = (b_i - sum_upper) / A_ii
        Solution_nPlus1.SetValue(target_j, target_i, (b_i - sum_upper) / A_ii);
    }
}


std::size_t GaussianElimination::FindPivotRowIndex(
    const CoefficientMatrix& A,
    const std::size_t PivotIndex)
{
    const std::size_t NumARows = A.GetNumRows();

    // Ensure the pivot index is within matrix bounds
    if (PivotIndex >= NumARows)
    {
        throw std::out_of_range("FindPivotRowIndex: PivotIndex is out of range.");
    }

    // Initialize the best candidate with the current diagonal element
    std::size_t CandidateRowIndex = PivotIndex;
    double MaxCandidateAbsValue = std::fabs(A.GetValue(PivotIndex, PivotIndex));

    // Search subsequent rows for the largest absolute pivot value
    for (std::size_t RowIndex = PivotIndex + 1; RowIndex < NumARows; ++RowIndex)
    {
        const double CurrentAbsValue = std::fabs(A.GetValue(RowIndex, PivotIndex));

        if (CurrentAbsValue > MaxCandidateAbsValue)
        {
            MaxCandidateAbsValue = CurrentAbsValue;
            CandidateRowIndex = RowIndex;
        }
    }

    return CandidateRowIndex;
}



void GaussianElimination::SwapRows(
    CoefficientMatrix& A,
    RHS& RHS,
    std::size_t FirstRowIndex,
    std::size_t SecondRowIndex)
{
    const std::size_t NumACols = A.GetNumColumns();

    if (FirstRowIndex >= A.GetNumRows() || SecondRowIndex >= A.GetNumRows())
    {
        throw std::out_of_range("SwapRows: row index is out of range.");
    }

    if (RHS.GetSize() != A.GetNumRows())
    {
        throw std::invalid_argument("SwapRows: size of RHS must match matrix size.");
    }

    if (FirstRowIndex == SecondRowIndex)
    {
        return;
    }

    for (std::size_t ColIndex = 0; ColIndex < NumACols; ++ColIndex)
    {
        const double Temp = A.GetValue(FirstRowIndex, ColIndex);
        A.SetValue(FirstRowIndex, ColIndex, A.GetValue(SecondRowIndex, ColIndex));
        A.SetValue(SecondRowIndex, ColIndex, Temp);
    }

    const double TempRHS = RHS.GetValue(FirstRowIndex);
    RHS.SetValue(FirstRowIndex, RHS.GetValue(SecondRowIndex));
    RHS.SetValue(SecondRowIndex, TempRHS);
}
