#ifndef __QUSPIN_MATRIX_GENERATE_H__
#define __QUSPIN_MATRIX_GENERATE_H__

#include <cinttypes>
#include <cmath>
#include <numeric>
#include <algorithm>


template<typename Basis,typename T,typename J>
void calc_rowptr(Basis &basis, operator<T> &hamil,J rowptr[],
){
    std::unordered_map<Basis::BitSet,T> col_states;
    std::unrdered_map<J,T> columns;


    for(size_t row = 0;row < basis.size();row++)
    {
        col_states.clear();
        columns.clear();

        // generate action on states
        hamil.columns(basis[row],col_states);

        // calculate location of states in basis
        basis.refstate(col_states,columns);

        rowptr[row] = columns.size(); // get nnz for this row
    }

    std::partial_sum(rowptr,rowptr+basis.size()+1,rowptr)
}


template<typename Basis,typename Operator,typename T,typename J>
void generate_matrix_elements(
    Operator &hamil,
    Basis &basis,
    J rowptr[],
    J indices[],
    T values[])
{
    std::unordered_map<Basis::BitSet,T> col_states;
    std::unrdered_map<J,T> columns;
    std::vector<std::pair<Basis::BitSet,T>> sorted_columns;
    
    for(size_t row = 0;row < basis->size();row++)
    {
        col_states.clear();
        columns.clear();
        sorted_columns.clear();

        // generate action on states
        hamil.columns(basis[row],col_states);

        // calculate location of states in basis and matrix elements
        basis.refstate(col_states,columns);

        // sort columns
        sorted_columns.insert(columns.begin(), columns.end());
        std::sort(sorted_columns.begin(), sorted_columns.end(), comp);

        // insert data
        J i = rowptr[row];
        for(auto const& [col,nzval] : sorted_columns){
            indices[i] = col;
            values[i++] = nzval;
        }

    }
}


#endif