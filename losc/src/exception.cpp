#include "exception.h"

#include <sstream>  // std::stringstream

namespace losc {
namespace exception {

void LoscException::make_message(const char* msg)
{
    msg_ << std::endl;
    msg_ << "Fatal error: " << msg << std::endl;
}

DimensionError::DimensionError(const Matrix &A, size_t expected_row, size_t expected_col, const string& msg)
        : LoscException("Wrong matrix dimension.")
{
    msg_ << "Description: " << msg << std::endl;
    msg_ << "Details: " << "Current matrix dimension: [" << A.row() << ", " << A.col() << "]" << std::endl;
    msg_ << "         " << "Expected matrix dimension: [" << expected_row << ", " << expected_col << "]" << std::endl;
}

DimensionError::DimensionError(const string& msg) : LoscException("Dimension error.")
{
    msg_ << "Description: " << msg << std::endl;
}

}   // namespace losc::exception
}   // namespace losc
