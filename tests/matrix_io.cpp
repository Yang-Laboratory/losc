#include <assert.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <unistd.h>

#include "matrix_io.hpp"
#include <losc/exception.hpp>

namespace test {
/**
 * @note The txt file `fname` will always be overwritten if `Mat` is not empty.
 */
void write_matrices_to_txt(vector<MatrixXd> &Mat, const string &fname,
                           size_t num_per_line)
{
    if (Mat.size() == 0)
        return;

    FILE *f = fopen(fname.c_str(), "w");
    if (f == NULL)
        throw losc::exception::LoscException(
            "Cannot open file to write matrices:" + fname);
    for (size_t i = 0; i < Mat.size(); ++i) {
        const MatrixXd &A = Mat[i];
        if (i == 0)
            fprintf(f, "Dimension,%zu,%zu\n", A.rows(), A.cols());
        else
            fprintf(f, "\nDimension,%zu,%zu\n", A.rows(), A.cols());

        // loop over each elements.
        size_t n = 0;
        const size_t size = A.size();
        for (size_t ii = 0; ii < A.rows(); ++ii) {
            for (size_t jj = 0; jj < A.cols(); ++jj) {
                const auto val = A(ii, jj);
                if (n != num_per_line && ii != size - 1) {
                    fprintf(f, "%.16e,", val);
                    ++n;
                } else if (n == num_per_line && ii != size - 1) {
                    fprintf(f, "%.16e\n", val);
                    n = 0;
                } else {
                    fprintf(f, "%.16e", val);
                }
            }
        }
    }
    fclose(f);
}

/**
 * @brief split string with given delimeter.
 * @param [in] str: input string to be splitted.
 * @param [in] delim: delimeter.
 * @param [in, out] str_split: On exit, it stores the splitted string.
 */
static void str_split(std::string &str, const char delim,
                      std::vector<std::string> &str_split)
{
    str_split.clear();
    std::stringstream str_stream(str + std::string(1, delim));
    std::string word;
    while (std::getline(str_stream, word, delim)) {
        str_split.push_back(word);
    }
}

std::vector<MatrixXd> read_matrices_from_txt(const string &fname)
{
    std::ifstream fin;
    fin.open(fname);
    if (!fin)
        throw losc::exception::LoscException(
            "Cannot open file to write matrices:" + fname);

    bool is_first_line = true;
    std::string line;
    std::vector<std::string> line_split;
    std::vector<MatrixXd> rst;
    MatrixXd current_matrix;
    std::size_t n_row = 0;
    std::size_t n_col = 0;
    std::size_t count = 0;

    // loop over each line in the txt file
    std::size_t line_num = 0;
    while (std::getline(fin, line)) {
        str_split(line, ',', line_split);
        auto p_data = line_split.begin();
        if (line_split[0] == "Dimension") {
            // now it starts to read a new matrix.
            // Before that, check if the last matrix is read successfully or
            // not.
            if (!is_first_line) {
                if (count != current_matrix.size()) {
                    throw losc::exception::LoscException(
                        "Error in read matrix, unmatched element size: file "
                        "name = " +
                        fname);
                }
            }
            // read dimension information for the new matrix.
            if (line_split.size() < 3) {
                throw losc::exception::LoscException(
                    "Cannot read matrix dimension: file name = " + fname);
            }
            try {
                n_row = std::stoi(line_split[1]);
                n_col = std::stoi(line_split[2]);
            } catch (...) {
                throw losc::exception::LoscException(
                    "Failed to read matrix dimension. file name = " + fname);
            }
            // create a new matrix and push it back to matrix vector.
            current_matrix = MatrixXd(n_row, n_col);
            current_matrix.setZero();
            rst.push_back(current_matrix);
            p_data += 3;
            count = 0; // set matrix element count as zero.
        }
        // put matrix element in current line into the matrix object.
        for (; p_data != line_split.end(); ++p_data) {
            try {
                const size_t i = count / n_col;
                const size_t j = count % n_col;
                // current_matrix->data()[count] = std::stod(*p_data);
                current_matrix(i, j) = std::stod(*p_data);
                ++count;
            } catch (...) {
                std::stringstream msg;
                msg << "Failed to read matrix element at line:" << line_num + 1
                    << ", column:" << p_data - line_split.begin() + 1
                    << std::endl;
                msg << "Failed element value: " << *p_data << std::endl;
                msg << "File name: " << fname << std::endl;
                throw losc::exception::LoscException(msg.str());
            }
        }
        // at this point, the first line has been parsed.
        is_first_line = false;
        line_num++;
    }
    fin.close();
    return rst;
}

} // namespace test
