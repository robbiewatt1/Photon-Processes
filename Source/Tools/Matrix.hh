#ifndef Matrix_hh
#define Matrix_hh

#include <iostream>
#include <cassert>
#include <iomanip>

#include "H5Cpp.h"


template <typename T>
class Matrix
{
    private:
        int m_nRow;
        int m_nColumn;
        T *m_data;
    
    public:

        Matrix():
        m_nRow(0), m_nColumn(0), m_data(nullptr)
        {
        }

        Matrix(int nRow, int nColumn):
        m_nRow(nRow), m_nColumn(nColumn)
        {
            // allocate the array space
            if (m_nRow * m_nColumn > 0)
            {
                m_data = new T [m_nRow * m_nColumn];
            } else
            {
                m_data = nullptr;
            }
        }

        Matrix(const std::initializer_list<std::initializer_list<T>> list):
        m_nRow(list.size()), m_nColumn(list.begin()->size())
        {
            m_data = new T [m_nRow * m_nColumn];
            int i = 0, j = 0;
            for (const auto& l : list)
            {
                for (const auto& v : l)
                {
                    m_data[i + j * m_nRow]
                }
            }
        }

        // Copy constructor 
        Matrix(const Matrix &matrix)
        {
            m_nRow = matrix.m_nRow;
            m_nColumn = matrix.m_nColumn;
            m_data = new T [matrix.m_nRow * matrix.m_nColumn];

            for (int i = 0; i < matrix.m_nRow * matrix.m_nColumn; i++)
            {
                m_data[i] = matrix.m_data[i];
            }
        }

        ~Matrix()
        {
            if (m_data)
            {
                delete [] m_data;
            }        
        }

        Matrix<T> deepCopy() const
        {
            Matrix<T> copy(m_nRow, m_nColumn);
            for (int i = 0; i < m_nRow * m_nColumn; i++)
            {

                copy.m_data[i] = m_data[i];
            }
            return copy;
        }
        
        int getRows() const
        {
            return m_nRow;
        }

        int getColumns() const
        {
            return m_nColumn;
        }

        T* begin() const
        {
            return m_data;
        }

        T* end() const
        {
            return &m_data[m_nRow * m_nColumn - 1];
        }

        // Prints Matrix to the screen
        void print() const
        {
            std::cout << *this;
        }

        void empty()
        {
            for (int i = 0; i < m_nRow * m_nColumn; i++)
            {
                m_data[i] = 0.0;
            }
        }
        
	void save(const H5std_string fileName, const H5std_string dataName)
        {
            H5::Exception::dontPrint(); 
            H5::H5File* file;
            try
            {
                file = new H5::H5File(fileName, H5F_ACC_RDWR);
            } catch (const H5::Exception& e)
            {
                file = new H5::H5File(fileName, H5F_ACC_EXCL);
            }

            int dataRank = 2;
            hsize_t dataDim[2];
            dataDim[0] = m_nRow;
            dataDim[1] = m_nColumn;
            H5::DataSpace dataSpace(dataRank, dataDim);
            try
            {
                H5::DataSet* dataSet = new H5::DataSet(
                    file->createDataSet(dataName, H5::PredType::NATIVE_DOUBLE,
                        dataSpace));
                dataSet->write(m_data, H5::PredType::NATIVE_DOUBLE);
                file->close();
                delete file;
                delete dataSet;
            } catch (const H5::Exception& e)
            {
                std::cerr << "Error: Data set: \'" << dataName
                    << "\'' probabily already exists." << std::endl;
                delete file;
                exit(1);
            }
        }

        void open(const H5std_string fileName, const H5std_string dataName)
        {
            // Clear current data
            if (m_data)
            {
                delete [] m_data;
            }

            H5::H5File* file = new H5::H5File(fileName, H5F_ACC_RDONLY);
            H5::DataSet dataset = file->openDataSet(dataName);
            H5::DataSpace dataspace = dataset.getSpace();
            int rank = dataspace.getSimpleExtentNdims();

            if(rank != 2)
            {
                std::cerr << "Error: Data set " << dataName << " has rank: " << rank
                    << ". This cannot be opend as a matrix." << std::endl;
                exit(1);
            }
            hsize_t dims[2];
            dataspace.getSimpleExtentDims(dims, NULL);
            H5::DataSpace mSpace(2, dims);

            // Allocate new space
            m_data = new T [dims[0] * dims[1]];
            m_nRow = dims[0];
            m_nColumn = dims[1];
            H5::DataType datatype = H5Dget_type(dataset. getId());
            dataset.read(m_data, datatype, mSpace, dataspace);
            delete file;
        }
        
	// Sums all the elements together
        T sum()
        {
            T sum(0);
            for (int i = 0; i < m_nRow * m_nColumn; i++)
            {

                sum += m_data[i];
            }
            return sum;
        }

        // Appends a row to the bottom of the matrix
        void appendRows(const Matrix<T> &row)
        {
            assert(m_nColumn == row.GetColumns());
            //  assert(row.GetColumns == m_nColumn);
            T *dataNew = new T [(m_nRow + row.GetRows()) * row.GetColumns()];

            // fill the new array
            for (int i = 0; i < m_nRow * m_nColumn; i++)
            {
                dataNew[i] = m_data[i];
            }

            for (int i = 0; i < row.GetRows() * row.GetColumns(); i++)
            {
                dataNew[i + m_nRow * m_nColumn] = row.Iter()[i];
            }

            delete [] m_data;
            m_data = dataNew;
            m_nRow = m_nRow + row.GetRows();
            m_nColumn = row.GetColumns();
        }

    public:

        // return the value of the matrix at point (rowIndex,columnIndex). Much prefered over 
        // [rowIndex,columnIndex] as preforms checks that index is in bound.
        T& operator()(int rowIndex, int columnIndex)
        {
            assert(rowIndex >= 0 && rowIndex < m_nRow);
            assert(columnIndex >= 0 && columnIndex < m_nColumn);

            return m_data[rowIndex * m_nColumn + columnIndex];
        }

        T& operator()(int rowIndex, int columnIndex) const
        {
            assert(rowIndex >= 0 && rowIndex < m_nRow);
            assert(columnIndex >= 0 && columnIndex < m_nColumn);

            return m_data[rowIndex * m_nColumn + columnIndex];
        }

        // return the column at the given row index. This method is no recomeneded as the are 
        // no proper checks to determine if the index is out of range of the array
        T* operator[](const int rowIndex)
        {
            return &m_data[rowIndex * m_nColumn];
        }

        T* operator[](int rowIndex) const
        {
            return &m_data[rowIndex * m_nColumn];
        }

        Matrix<T>& operator=(const Matrix<T> &matrix)
        {
            // Self assigmnet gard
            if( this == &matrix)
            {
                return *this;
            }

            // Delete any data in new matrix is holding
            if (m_data)
            {
                delete[] m_data;
            }

            // Copy variables over
            m_nRow = matrix.m_nRow;
            m_nColumn = matrix.m_nColumn;
            if(matrix.m_data)
            {
                // allocate the array space
                m_data = new T [m_nRow * m_nColumn];
                // Copy the data over
                for(int i = 0; i < m_nRow * m_nColumn; i++)
                {

                    m_data[i] = matrix.m_data[i];
                }
            }else
            {
                m_data = 0;
            }

            return *this;
            // Copys a matirx to a new variable
        }

    friend Matrix<T> operator*(const Matrix<T> &matrix1, const Matrix<T> &matrix2)
    {
        Matrix<T> newMatrix = Matrix<T>(matrix1.m_nRow,matrix2.m_nColumn);

        for(int i = 0; i < matrix1.m_nRow; i++)
        {
            for(int j = 0; j < matrix2.m_nColumn; j++)
            {
                double elementValue = 0.0;
                for(int k = 0; k < matrix1.m_nColumn; k++)
                {
                    elementValue = elementValue + matrix1(i, k) * matrix2(k, j);
                }
                newMatrix(i, j) = elementValue;
            }
        }
        return newMatrix;
    }

    friend Matrix<T> operator*(const T &scalar, const Matrix<T> &matrix)
    {
        Matrix<T> newMatrix = Matrix<T>(matrix.m_nRow,matrix.m_nColumn);

        for(int i = 0; i < matrix.m_nRow; i++)
        {
            for(int j = 0; j < matrix.m_nColumn; j++)
            {
                newMatrix(i,j) = scalar * matrix(i, j);
            }
        }
        return newMatrix;
    }

    friend Matrix<T> operator+(const Matrix<T> &matrix1, const Matrix<T> &matrix2)
    {
        assert(matrix1.m_nRow == matrix2.m_nRow);
        assert(matrix1.m_nColumn == matrix2.m_nColumn);

        Matrix<T> newMatrix = Matrix<T>(matrix1.m_nRow,matrix1.m_nColumn);
        for(int i = 0; i < matrix1.m_nRow; i++)
        {
            for(int j = 0; j < matrix1.m_nColumn; j++)
            {
            newMatrix(i,j) = matrix1(i, j) + matrix2(i, j);
            }
        }
        return newMatrix;
    }

    friend Matrix<T> operator-(const Matrix<T> &matrix1, const Matrix<T> &matrix2)
    {
        assert(matrix1.m_nRow == matrix2.m_nRow);
        assert(matrix1.m_nColumn == matrix2.m_nColumn);

        Matrix<T> newMatrix = Matrix<T>(matrix1.m_nRow,matrix1.m_nColumn);
        for(int i = 0; i < matrix1.m_nRow; i++)
        {
            for(int j = 0; j < matrix1.m_nColumn; j++)
            {
            newMatrix(i,j) = matrix1(i, j) - matrix2(i, j);
            }
        }
        return newMatrix;
    }

    friend std::ostream& operator << (std::ostream& out, const Matrix<T>& matrix)
    {
        out << "[ ";
        for(int i = 0; i <= matrix.m_nRow * matrix.m_nColumn - 1; i++)
        {
            out << std::fixed << matrix.m_data[i] << ", ";
            if((i + 1) % matrix.m_nColumn == 0 && i != matrix.m_nRow * matrix.m_nColumn - 1)
            {
                out << "\n  ";    
            }else if(i == matrix.m_nRow * matrix.m_nColumn - 1)
            {
                out << "]" << std::endl;
            }
        }
        return out;
    }
};
#endif
