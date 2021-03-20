#include <losc/exception.hpp>
#include <sstream> // std::stringstream

namespace losc {
namespace exception {

void LoscException::make_message(const char *msg)
{
    stringstream s;
    s << std::endl;
    s << "Fatal error: " << msg << std::endl;
    msg_ = s.str();
}

DimensionError::DimensionError(ConstRefMat &A, size_t expected_row,
                               size_t expected_col, const string &msg)
    : LoscException("Wrong matrix dimension.")
{
    stringstream s;
    s << "Description: " << msg << std::endl;
    s << "Details: "
         << "Current matrix dimension: [" << A.rows() << ", " << A.cols() << "]"
         << std::endl;
    s << "         "
         << "Expected matrix dimension: [" << expected_row << ", "
         << expected_col << "]" << std::endl;
    msg_ = s.str();
}

DimensionError::DimensionError(const string &msg)
    : LoscException("Dimension error.")
{
    stringstream s;
    s << "Description: " << msg << std::endl;
    msg_ = s.str();
}

} // namespace exception
} // namespace losc
