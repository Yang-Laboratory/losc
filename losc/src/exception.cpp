/**
 * @file
 * @brief Losc library exception realted definition.
 */
#include "exception.h"
#include <sstream> // std::stringstream

namespace losc {
namespace exception {

void LoscException::make_message(const char *msg)
{
    msg_ << std::endl;
    msg_ << "Fatal error: " << msg << std::endl;
}

DimensionError::DimensionError(ConstRefMat &A, size_t expected_row,
                               size_t expected_col, const string &msg)
    : LoscException("Wrong matrix dimension.")
{
    msg_ << "Description: " << msg << std::endl;
    msg_ << "Details: "
         << "Current matrix dimension: [" << A.rows() << ", " << A.cols() << "]"
         << std::endl;
    msg_ << "         "
         << "Expected matrix dimension: [" << expected_row << ", "
         << expected_col << "]" << std::endl;
}

DimensionError::DimensionError(const string &msg)
    : LoscException("Dimension error.")
{
    msg_ << "Description: " << msg << std::endl;
}

} // namespace exception
} // namespace losc
