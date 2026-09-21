#include "Classes2D/FileWriter2D/FileWriter2D.hpp"

void FileWriter2D::WriteField2D(
    Field2D& Field2D_Obj,  
    Mesh2D& Mesh2D_Obj ,
    FieldType FieldType_,
    std::string RelativePath)
{
    std::string FileName;
    if (FieldType_ == FieldType::Initial_Numerical)
    {
        std::cout<<"Initialization Completed, Writing The Initial Field2D To File..."<<std::endl;
        FileName = "Initial_Numerical2D";
    }
    else if (FieldType_ == FieldType::Steady_State_Numerical)
    {
        std::cout<<"Writing The Steady State Numerical Field2D To File..."<<std::endl;
        FileName = "Steady_State_Numerical2D";
    }

    std::string FileExtension = "dat";
    std::filesystem::path dirPath(RelativePath);
    if (!std::filesystem::exists(dirPath))
    {
        std::filesystem::create_directories(dirPath);
    }

    std::filesystem::path fullFilePath = dirPath / (FileName + "." + FileExtension);

    std::ofstream file(fullFilePath);


    if (file.is_open()) 
    {
    const std::size_t nrows = Field2D_Obj.GetNrows();
    const std::size_t ncols = Field2D_Obj.GetNcols();
        file << std::fixed << std::setprecision(6);
    for (std::size_t row = 0; row < nrows; ++row)
        {    
        //const std::size_t yIndex = nrows - 1 - row;
        for (std::size_t col = 0; col < ncols; ++col)
            {    
                file << std::left
                 << std::setw(20) << Mesh2D_Obj.GetXPositions(col)
                 << std::setw(20) << Mesh2D_Obj.GetYPositions(row)
                 << std::setw(20) << Field2D_Obj.GetValue(row, col)
                 << '\n';
            }
        }
                             
    std::cout << FileName << "." << FileExtension << " Created Successfully At: " << fullFilePath.string() << std::endl;
    file.close();
    }

    else 
    {
        std::cout << FileName << "." << FileExtension << " Failed To Create Because The Path Folder Is WRONG." << std::endl;
    } 

}

void FileWriter2D::WriteField2D(
    AnalyticalDiffusion2D& Analytical2D_Obj,  
    Mesh2D& Mesh2D_Obj ,
    [[maybe_unused]] FieldType FieldType_,
    std::string RelativePath)
{
    std::cout<<"Analytical Tempreture values have been calculated, Writing The Analytical velocity File"<<std::endl;
    std::string FileName = "Steady_State_Analytical2D";
    std::string FileExtension = "dat";
    std::filesystem::path dirPath(RelativePath);
    if (!std::filesystem::exists(dirPath))
    {
        std::filesystem::create_directories(dirPath);
    }

    std::filesystem::path fullFilePath = dirPath / (FileName + "." + FileExtension);

    std::ofstream file(fullFilePath);

    if (file.is_open()) 
    {
    const std::size_t nrows = Analytical2D_Obj.GetNrows();
    const std::size_t ncols = Analytical2D_Obj.GetNcols();
        file << std::fixed << std::setprecision(6);
    for (std::size_t row = 0; row < nrows; ++row)
        {    
        //const std::size_t yIndex = nrows - 1 - row;
        for (std::size_t col = 0; col < ncols; ++col)
            {    
                file << std::left
                 << std::setw(20) << Mesh2D_Obj.GetXPositions(col)
                 << std::setw(20) << Mesh2D_Obj.GetYPositions(row)
                 << std::setw(20) << Analytical2D_Obj.GetValue(row, col)
                 << '\n';
            }
        }
                             
    std::cout << FileName << "." << FileExtension << " Created Successfully At: " << fullFilePath.string() << std::endl;
    file.close();
    }

    else 
    {
        std::cout << FileName << "." << FileExtension << " Failed To Create Because The Path Folder Is WRONG." << std::endl;
    } 

}


void FileWriter2D::WriteFullSystemCSV(
    const CoefficientMatrix& A,
    const RHS& RHSD1D,
    const Mesh2D& Mesh2D_Obj,
    std::string RelativePath)
{
    namespace fs = std::filesystem;

    // Ensure the output directory exists
    if (!RelativePath.empty() && !fs::exists(RelativePath))
    {
        fs::create_directories(RelativePath);
    }

    std::string fullFilePath = RelativePath + "/Full_Linear_System.csv";
    std::ofstream outFile(fullFilePath);

    if (!outFile.is_open())
    {
        throw std::runtime_error("FileWriter2D::WriteFullSystemCSV: Failed to open file " + fullFilePath);
    }

    int nrows = Mesh2D_Obj.GetNumOfYNodes();
    int ncols = Mesh2D_Obj.GetNumOfXNodes();
    int totalNodes = nrows * ncols;

    // Write CSV Header
    outFile << "Row_Index,Node(j;i),Boundary_Type,RHS_Value";
    for (int col = 0; col < totalNodes; ++col)
    {
        int c_j = col / ncols;
        int c_i = col % ncols;
        outFile << ",Col_" << col << " (" << c_j << ";" << c_i << ")";
    }
    outFile << "\n";

    // Set floating-point precision for accurate inspection
    outFile << std::scientific << std::setprecision(6);

    // Iterate through every global row
    for (int row = 0; row < totalNodes; ++row)
    {
        int j = row / ncols;
        int i = row % ncols;

        // Determine node type following the strict priority order: Top -> Bottom -> Right -> Left -> Interior
        std::string boundaryType = "Interior";
        if (j == nrows - 1)
        {
            boundaryType = "Top";
        }
        else if (j == 0)
        {
            boundaryType = "Bottom";
        }
        else if (i == ncols - 1)
        {
            boundaryType = "Right";
        }
        else if (i == 0)
        {
            boundaryType = "Left";
        }

        // Write row metadata and RHS value
        outFile << row << ",\"" << "(" << j << "," << i << ")\"," << boundaryType << "," << RHSD1D.GetValue(row);

        // Write each matrix coefficient in this row
        for (int col = 0; col < totalNodes; ++col)
        {
            outFile << "," << A.GetValue(row, col);
        }
        outFile << "\n";
    }

    outFile.close();
    std::cout << "Successfully written full system to: " << fullFilePath << std::endl;
}
