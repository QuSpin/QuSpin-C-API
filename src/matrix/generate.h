#ifndef __QUSPIN_MATRIX_GENERATE_H__
#define __QUSPIN_MATRIX_GENERATE_H__

#include <cinttypes>
#include <cmath>

template<typename Basis,typename T,typename J>
void calc_rowptr(Basis &basis, operator<T> &hamil,J rowptr[],
){
    J nnz = 0;
    for(size_t row = 0;row < basis.size();row++)
    {
        rowptr[row] = nnz;
        nnz += hamil.nnz_columns(basis,row); // get nnz for this row
    }
    rowptr[basis.size()] = nnz;
}


template<typename Basis,typename T,typename J>
void generate_matrix(
    Basis &basis,
    operator<T> &hamil,
    J rowptr[],
    J indices[],
    T values[])
{
    for(size_t row = 0;row < basis->size();row++)
    {
        auto results = hamil.columns(basis,row);
        auto results_iter = results.begin();
        for(J ind=rowptr[row],cind=0;ind<rowptr[row+1];++ind,++cind){
            indices[ind] = results_iter.first();
            values[ind] = results_iter.second();
        }
    }
}


#endif