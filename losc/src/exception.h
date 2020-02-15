#ifndef _LOSC_EXCEPTION_H_
#define _LOSC_EXCEPTION_H_

#include <stdexcept>
#include <string>
#include <sstream>  // std::stringstream
#include <matrix/matrix.h>

namespace losc {

namespace exception {

using matrix::Matrix;
using std::string;
using std::stringstream;

class LoscException : public std::runtime_error
{
    private:
    void make_message(const char* msg);

    protected:
    stringstream msg_;

    public:
    LoscException(const std::string& msg) : std::runtime_error(msg) {make_message(msg.c_str());}

    const char* what() const noexcept override {return msg_.str().c_str();}
};

class DimensionError : public LoscException
{
    public:
    DimensionError(const Matrix &A, size_t expected_row, size_t expected_col, const string& msg);

    DimensionError(const string& msg);
};

}
}


#endif // _LOSC_EXCEPTION_H_
