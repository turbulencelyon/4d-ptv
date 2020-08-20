/*
 * Simple CMatrix object.
 *
 * Apparently taken from here: http://www.softwareandfinance.com/CPP/Matrix_Determinant.html
 *
 */

#ifndef STM_CMATRIX_H
#define STM_CMATRIX_H

class CMatrix {

    private:
        int m_rows;
        int m_cols;
        CMatrix();
    
    public:
        double **m_pData;

    CMatrix(int rows, int cols);
    CMatrix(const CMatrix &other);
    ~CMatrix();
    
    void GetInput();
    void FillSimulatedInput();
    double Determinant();
    CMatrix& operator =(const CMatrix &other);
    CMatrix CoFactor();
    CMatrix Adjoint();
    CMatrix Transpose();
    CMatrix Inverse();
    CMatrix operator +(const CMatrix &other);
    CMatrix operator -(const CMatrix &other);
    CMatrix operator *(const CMatrix &other);
    bool operator ==(const CMatrix &other);
    friend std::istream& operator >>(std::istream &is, CMatrix &m);
    friend std::ostream& operator <<(std::ostream &os, const CMatrix &m);
};

std::istream& operator >>(std::istream &is, CMatrix &m);
std::ostream& operator <<(std::ostream &os, const CMatrix &m);

#endif
