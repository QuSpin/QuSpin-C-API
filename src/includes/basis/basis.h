#ifndef __QUSPIN_BASIS_BASIS_H__
#define __QUSPIN_BASIS_BASIS_H__

#include <unordered_map>
#include <algorithm>

#include "matrix/operater.h"

namespace quspin::basis {

template<typename space_t, typename symmetry_t>
class basis
{
private:
    symmetry_t symmetry; // OK to copy
    space_t * space; // ? on whether to copy

public:

    typedef space_t::bitset_t bitset_t;
    typedef space_t::index_t index_t;
    typedef space_t::norm_t norm_t;

    basis(symmetry_t& _symmetry, space_t& _sapce) : symmetry(_symmetry), space(&_space) {}
    ~basis() {}

    template<typename T>
    void ref_states_conj(
        const typename index_t i,
        std::unrdered_map<typename bitset_t,T> col_states, 
        std::unrdered_map<typename index_t,T>& columns
    ){
        for(auto const& ele : col_states){
            auto refstate = symmetry.get_refstate(ele.first);
            typename index_t j = state.get_index(refstate.first);
            typename norm_t norm_j = space.get_norm(j);
            typename norm_t norm_i = space.get_norm(i);
            
            columns[state_index] = ele.second * std::conj(refstate.second) * std::sqrt(double(norm_j)/ norm_i);
        }
    }

    template<typename T>
    void ref_states(
        const typename index_t i,
        std::unrdered_map<typename bitset_t,T> row_states, 
        std::unrdered_map<typename index_t,T>& rows    
    ){
        for(auto const& ele : row_states){
            auto refstate = symmetry.get_refstate(ele.first);
            typename index_t j = state.get_index(refstate.first);
            typename norm_t norm_j = space.get_norm(j);
            typename norm_t norm_i = space.get_norm(i);
            
            rows[state_index] = ele.second * refstate.second * std::sqrt(double(norm_j)/ norm_i);
        }
    }

    template<typename T>
    void calc_rowptr(
        std::vector<operator_string<T>>& pterms,
        std::vector<dense_term<T>>& dterms,
        typename index_t rowptr[]
    ){

        J n_row = space -> size();

        std::unordered_map<typename bitset_t,T> col_states;
        std::unordered_map<typename index_t,T> columns;

        rowptr[0] = 0;
        for(typename index_t row = 0;row < n_row;++row)
        {
            col_states.clear();
            columns.clear();
            auto state = space.get_state(row);
            // generate action on states
            for(auto const& pterm : pterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            for(auto const& dterm : dterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            // calculate location of states in basis
            this->ref_state_conj(col_states,columns);
            // insert number of non-zeros elements for this row
            rowptr[row] = columns.size();
        }

        typename index_t nnz = 0;
        for(typename index_t row = 0;row < n_row;++row){
            typename index_t tmp = rowptr[row];
            rowptr[row] = nnz;
            nnz += tmp;
        }
        rowptr[n_row+1] = nnz;

    }

    template<typename T>
    void calc_rowptr(
        std::vector<operator_string<T>>& pterms,
        std::vector<dense_term<T>>& dterms,
        typename index_t rowptr[],
        typename index_t indices[],
        T values[]
    ){

        J n_row = space -> size();

        std::unordered_map<typename bitset_t,T> col_states;
        std::unordered_map<typename index_t,T> columns;
        std::vector<std::pair<typename index_t,T>> sorted_columns;


        rowptr[0] = 0;
        for(typename index_t row = 0;row < n_row;++row)
        {
            col_states.clear();
            columns.clear();
            sorted_columns.clear();

            auto state = space.get_state(row);
            // generate action on states
            for(auto const& pterm : pterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            for(auto const& dterm : dterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            // calculate location of states in basis
            this->ref_state_conj(col_states,columns);

            // sort columns
            sorted_columns.insert(columns.begin(), columns.end());
            std::sort(sorted_columns.begin(), sorted_columns.end(), 
                [](std::pair<typename index_t,T> lhs,std::pair<typename index_t,T> rhs) -> bool 
                {
                    return lhs.first < rhs.first;
                }
            );

            // insert data
            typename index_t i = rowptr[row];
            for(auto const& [col,nzval] : sorted_columns){
                indices[i] = col;
                values[i++] = nzval;
            }

        }

    }


};

}

#endif