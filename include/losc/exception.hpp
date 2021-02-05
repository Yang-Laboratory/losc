/**
 * @file exception.hpp
 * @brief C++ interface for the LOSC exceptions.
 */

#ifndef _LOSC_INCLUDE_LOSC_EXCEPTION_HPP_
#define _LOSC_INCLUDE_LOSC_EXCEPTION_HPP_

#include <losc/eigen_def.hpp>
#include <sstream> // std::stringstream
#include <stdexcept>
#include <string>

namespace losc {

/**
 * @brief namespace for losc exceptions.
 */
namespace exception {

using losc::ConstRefMat;
using losc::ConstRefVec;
using losc::RefMat;
using losc::RefVec;
using std::string;
using std::stringstream;

class LoscException : public std::runtime_error {
  private:
    void make_message(const char *msg);

  protected:
    stringstream msg_;

  public:
    LoscException(const std::string &msg) : std::runtime_error(msg)
    {
        make_message(msg.c_str());
    }

    const char *what() const noexcept override { return msg_.str().c_str(); }
};

class DimensionError : public LoscException {
  public:
    DimensionError(ConstRefMat &A, size_t expected_row, size_t expected_col,
                   const string &msg);

    DimensionError(const string &msg);
};

} // namespace exception
} // namespace losc
#endif // _LOSC_SRC_EXCEPTION_H_
