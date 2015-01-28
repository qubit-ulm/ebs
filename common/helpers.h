#ifndef HELPERS_H
#define HELPERS_H

#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/type_traits/is_complex.hpp>

namespace helpers
{
    template<typename VectorType>
    inline void circshift(VectorType& v, int shift)
    {
        if (abs(shift) > v.n_elem)
            throw new std::invalid_argument("number of elements to shift must be less then vector length");

        VectorType tmp(abs(shift));

        if (shift >= 1)
            shiftLeft(v, tmp, shift);
        else
            shiftRight(v, tmp, abs(shift));
    }

    template<typename VectorType>
    inline void shiftLeft(VectorType& v, VectorType& tmp, int shift)
    {
        for (int i = 0; i < shift; i++)
            tmp(i) = v(i);

        for (int i = 0; i < v.n_elem - shift; i++)
            v(i) = v(i + shift);

        for (int i = 0; i < shift; i++)
            v(v.n_elem - shift + i) = tmp(i);
    }

    template<typename VectorType>
    inline void shiftRight(VectorType& v, VectorType& tmp, int shift)
    {
        for (int i = 0; i < shift; i++)
            tmp(shift - 1 - i) = v(v.n_elem - 1 - i);

        for (int i = v.n_elem - 1; i >= shift; i--)
            v(i) = v(i-shift);

        for (int i = 0; i < shift; i++)
            v(i) = tmp(i);
    }
}

namespace helpers
{
    namespace
    {
        std::tuple<int, int> getDimensions(std::istream& input)
        {
            using namespace std;
            string line;

            do {
                getline(input, line);
                assert(input.good());
            } while (boost::starts_with(line, "%"));

            int n = 0, col = 0;
            istringstream line_with_dims(line);

            line_with_dims >> n >> col;
            assert(n > 0 && col > 0);

            return make_tuple(n, col);
        }

        template <typename ScalarType>
        inline void readValue(const std::string &line, ScalarType &val)
        {
            using namespace std;

            istringstream newline(line);
            newline >> val;
        }

        template <typename ScalarType>
        inline void readValue(const std::string &line, std::complex<ScalarType> &val)
        {
            using namespace std;

            ScalarType re, img;
            istringstream newline(line);
            newline >> re >> img;
            val = complex<ScalarType>(re, img);
        }

        template <typename VectorType>
        int loadDataToVec(VectorType& vec, std::istream& input, const int n)
        {
            using namespace std;
            typedef typename VectorType::value_type ScalarType;

            int i = 0;
            string line;
            ScalarType val;

            while (getline(input, line) && i < n) {
                readValue(line, val);
                vec(i++) = val;
            }

            return i;
        }

        template <typename DataType>
        bool writeMMHeader(DataType& vec, std::ostream& output)
        {
            using namespace std;
            typedef typename DataType::value_type ScalarType;

            output.flags(ios_base::scientific);
            output.precision(16);

            if (boost::is_complex<ScalarType>::value)
                output << "%%MatrixMarket matrix array complex general\n";
            else
                output << "%%MatrixMarket matrix array real general\n";

            return true;
        }

        template <typename ScalarType>
        void writeValue(std::ostream& output, ScalarType value)
        {
            output << value << "\n";
        }

        template <typename ScalarType>
        void writeValue(std::ostream& output, std::complex<ScalarType> value)
        {
            output << value.real() << " " << value.imag() << "\n";
        }

        template <typename VectorType>
        int writeDataFromVec(const VectorType& vec, std::ostream &output)
        {
            using namespace std;
            int i;
            output << vec.size() << " " << 1 << "\n";

            for (i = 0; i < (int)vec.size(); i++)
                writeValue(output, vec(i));

            return i;
        }

	template <typename MatType>
	int writeDataFromMat(const MatType &mat, std::ostream &output)
	{
	   using namespace std;

	   size_t rows = mat.size1();
	   size_t cols = mat.size2();
	   output << rows  << " " << cols << " " << (cols * rows) << "\n";

	   size_t r, c;
	   for (r = 0; r < rows; r++) {
		   for (c = 0; c < cols; c++) {
		   	output << (r+1) << " " << (c+1) << " " << mat(r, c) << "\n";
		   }
	   }

	   return r * c;
	}
    }

    template <typename VectorType>
    bool loadMMVector(VectorType& vec, std::string filename)
    {
        using namespace std;

        ifstream input(filename.c_str(), ios::in);
        if (! input)
            return false;

        int n, col;
        std::tie(n, col) = getDimensions(input);

        vec.resize(n);
        int i = loadDataToVec(vec, input, n);

        input.close();

        if (i != n) {
            cerr << "Error while reading elements from file '" << filename << "'";
            return false;
        }

        return true;
    }

    template <typename VectorType>
    bool saveMMVector(const VectorType &vec, std::ostream &os)
    {
        writeMMHeader(vec, os);
        writeDataFromVec(vec, os);

        return true;
    }

    template <typename VectorType>
    bool saveMMVector(const VectorType &vec, std::string filename)
    {
        using namespace std;

        ofstream output(filename, ios::out);
        if (! output)
            return false;

        if (! saveMMVector(vec, output))
            return false;

        output.close();
        return true;
    }

    template <typename MatType>
    bool saveMMMatrix(const MatType &mat, std::ostream &os)
    {
	writeMMHeader(mat, os);
	writeDataFromMat(mat, os);

	return true;
    }

    template <typename MatType>
    bool saveMMMatrix(const MatType &mat, std::string filename)
    {
	using namespace std;

	ofstream output(filename, ios::out);
	if (! output)
	    return false;

	if (! saveMMMatrix(mat, output))
	    return false;

	output.close();
	return true;
    }
}

#endif // HELPERS_H
