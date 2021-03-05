/**
 * @file
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

/**
 * @brief The base exception in LOSC
 */
class LoscException : public std::runtime_error {
  private:
    void make_message(const char *msg);

  protected:
    stringstream msg_;

  public:
    /**
     * @brief Constructor of LoscException.
     * @param [in] msg the description of error.
     */
    LoscException(const std::string &msg) : std::runtime_error(msg)
    {
        make_message(msg.c_str());
    }

    /**
     * @brief Deconstructor of LoscException.
     */
    ~LoscException() {}

    const char *what() const noexcept override { return msg_.str().c_str(); }
};

/**
 * @brief The exception related to dimension error in LOSC
 */
class DimensionError : public LoscException {
  public:
    /**
     * @brief Constructor of DimensionError related to a matrix.
     * @param [in] A a LOSC matrix.
     * @param [in] expected_row the expected row number of `A`.
     * @param [in] expected_col the expected column number of `A`.
     * @param [in] msg the description of error.
     */
    DimensionError(ConstRefMat &A, size_t expected_row, size_t expected_col,
                   const string &msg);
    /**
     * @brief Constructor of DimensionError with a general message.
     * @param [in] msg the description of error.
     */
    DimensionError(const string &msg);

    /**
     * @brief Deconstructor of DimensionError.
     */
    ~DimensionError() {}
};

} // namespace exception
} // namespace losc
#endif // _LOSC_SRC_EXCEPTION_H_
