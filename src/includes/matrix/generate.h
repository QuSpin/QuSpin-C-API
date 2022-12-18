#ifndef __QUSPIN_MATRIX_GENERATE_H__
#define __QUSPIN_MATRIX_GENERATE_H__

#include <cinttypes>
#include <cmath>
#include <numeric>
#include <algorithm>

namespace quspin {

template<typename basis_t,typename T,typename J>
void calc_rowptr(basis_t &basis, operator<T> &hamil,J rowptr[],)
{
    J n_row = basis.size();
    std::unordered_map<typename basis_t::bitset_t,T> col_states;
    std::unrdered_map<J,T> columns;

    rowptr[0] = 0;
    for(J row = 0;row < n_row;row++)
    {
        col_states.clear();
        columns.clear();

        // generate action on states
        hamil.columns(basis[row],col_states);

        // calculate location of states in basis
        basis.ref_state_conj(col_states,columns);

        rowptr[row+1] = columns.size(); // get nnz for this row
    }

    std::partial_sum(rowptr,rowptr+n_row+1,rowptr)
}

template<typename basis_t,typename T>
void generate_matrix_elements(
    operator<T> &hamil,
    basis_t &basis,
    J rowptr[],
    J indices[],
    T values[])
{
    std::unordered_map<basis_t::bitset_t,T> col_states;
    std::unrdered_map<J,T> columns;
    std::vector<std::pair<typename basis_t::bitset_t,T>> sorted_columns;
    
    for(size_t row = 0;row < basis.size();row++)
    {
        col_states.clear();
        columns.clear();
        sorted_columns.clear();

        // generate action on states
        hamil.columns(basis[row],col_states);

        // calculate location of states in basis and matrix elements
        basis.ref_state_conj(row,col_states,columns);

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

}

#endif