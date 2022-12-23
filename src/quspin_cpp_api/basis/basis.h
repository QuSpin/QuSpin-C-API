#ifndef __QUSPIN_BASIS_BASIS_H__
#define __QUSPIN_BASIS_BASIS_H__

#include <unordered_map>
#include <algorithm>

#include "operater.h"

namespace quspin::basis {

template<typename subspace_t, typename symmetry_t>
class symmetric_basis
{
private:
    symmetry_t symmetry; // OK to copy
    shared_ptr<subspace_t> space;

public:

    typedef subspace_t::bitset_t bitset_t;
    typedef subspace_t::index_t index_t;
    typedef subspace_t::norm_t norm_t;

    basis(symmetry_t& _symmetry, shared_ptr<subspace_t> _sapce) : symmetry(_symmetry), space(_space) {}
    ~basis() {}

    shared_ptr<subspace_t> get_space() {return space;}
    symmetry_t get_symmetries() const {return symmetry;}

    template<typename J,typename T>
    void ref_states_conj(
        const J i,
        const std::unrdered_map<J,T>& col_states, 
        std::unrdered_map<J,T>& columns
    ) const {
        for(const auto& [state,raw_mat_ele] : col_states){
            const auto& [ref_state,charater] = symmetry.get_refstate(state);
            const J j = space->get_index(ref_state);
            const typename norm_t norm_j = space->get_norm(j);
            typename norm_t norm_i = space->get_norm(i);
            const T mat_ele =  raw_mat_ele * std::conj(charater) * std::sqrt(double(norm_j) / norm_i);
            (columns.contains(state_index) ?  columns[state_index] = mat_ele : columns[state_index] += mat_ele);
        }
    }

    template<typename J, typename T>
    void ref_states(
        const J i,
        const std::unrdered_map<J,T>& row_states, 
        std::unrdered_map<J,T>& rows    
    ) const {
        for(const auto& [state,raw_mat_ele] : row_states){
            const auto& [ref_state,charater] = symmetry.get_refstate(state);
            J j = space->get_index(ref_state);
            const auto norm_j = space->get_norm(j);
            const auto norm_i = space->get_norm(i);
            const T mat_ele =  raw_mat_ele * charater * std::sqrt(double(norm_j) / norm_i);
            (rows.contains(state_index) ? : rows[state_index] = mat_ele : rows[state_index] += mat_ele);
        }
    }

    template<typename J,typename T>
    void calc_rowptr(
        std::vector<operator_string<T>>& pterms,
        std::vector<dense_term<T>>& dterms,
        J rowptr[]
    ) const {

        J n_row = space -> size();

        std::unordered_map<typename bitset_t,T> col_states;
        std::unordered_map<J,T> columns;

        rowptr[0] = 0;
        for(J row = 0;row < n_row;++row)
        {
            col_states.clear();
            columns.clear();
            auto state = space->get_state(row);
            // generate action on states
            for(const auto& pterm : pterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            for(const auto& dterm : dterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            // calculate location of states in basis
            this->ref_state_conj<T>(col_states,columns);
            // insert number of non-zeros elements for this row
            rowptr[row] = columns.size();
        }

        J nnz = 0;
        for(J row = 0;row < n_row;++row){
            J tmp = rowptr[row];
            rowptr[row] = nnz;
            nnz += tmp;
        }
        rowptr[n_row+1] = nnz;

    }

    template<typename J,typename T>
    void calc_matrix(
        std::vector<operator_string<T>>& pterms,
        std::vector<dense_term<T>>& dterms,
        J rowptr[],
        J indices[],
        T values[]
    ) const {

        J n_row = space -> size();

        std::unordered_map<typename bitset_t,T> col_states;
        std::unordered_map<J,T> columns;
        std::vector<std::pair<J,T>> sorted_columns;


        rowptr[0] = 0;
        for(J row = 0;row < n_row;++row)
        {
            col_states.clear();
            columns.clear();
            sorted_columns.clear();

            auto state = space.get_state(row);
            // generate action on states
            for(const auto& pterm : pterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            for(const auto& dterm : dterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            // calculate location of states in basis
            this->ref_state_conj<T>(col_states,columns);

            // sort columns
            sorted_columns.insert(columns.begin(), columns.end());
            std::sort(sorted_columns.begin(), sorted_columns.end(), 
                [](std::pair<J,T> lhs,std::pair<J,T> rhs) -> bool 
                {
                    return lhs.first < rhs.first;
                }
            );

            // insert data
            J i = rowptr[row];
            for(const auto& [col,nzval] : sorted_columns){
                indices[i] = col;
                values[i++] = nzval;
            }

        }

    }

    template<typename T, typename X, typename Y>
    void on_the_fly(
        std::vector<operator_string<T>>& pterms,
        std::vector<dense_term<T>>& dterms,
        const Y a,
        const X * x, 
        const Y b, 
        Y  * y
    )const {

        if(b == Y(0.0)){
            std::fill(y,y+space->size(),0)
        }
        else{
            for(size_t index=0;index < space-<size();++index){ y[index] *= b; }
        }

        std::unordered_map<typename bitset_t,T> row_states;
        std::unordered_map<size_t,T> matrix_ele;

        for(J row=0;row < space->size();++row){
            row_states.clear();
            matrix_ele.clear();

            auto state = space.get_state(row);
            // generate action on states
            for(const auto& pterm : pterms){
                pterm.op<typename bitset_t>(state,col_states);
            }
            for(const auto& dterm : dterms){
                pterm.op<typename bitset_t>(state,col_states);
            }
            // calculate location of states in basis
            this->ref_state<T>(row_states,matrix_ele);

            Y total = 0;
            for(const auto& [col,nzval] : matrix_ele){total += nzval * x[col];}
            y[row] += a * total;

            matrix_ele.clear();
        }
    }

};

template<typename subspace_t>
class basis
{
private:
    shared_ptr<subspace_t> space;

public:

    typedef subspace_t::bitset_t bitset_t;
    typedef subspace_t::index_t index_t;
    typedef subspace_t::norm_t norm_t;

    basis(shared_ptr<subspace_t> _sapce) : space(_space) {}
    ~basis() {}

    shared_ptr<subspace_t> get_space() {return space;}


    template<typename J,typename T>
    void ref_states_conj(
        const J i,
        const std::unrdered_map<J,T>& col_states, 
        std::unrdered_map<J,T>& columns
    ) const {
        for(const auto& [state,raw_mat_ele] : col_states){
            const auto& [ref_state,charater] = symmetry.get_refstate(state);
            const J j = space->get_index(ref_state);
            const typename norm_t norm_j = space->get_norm(j);
            typename norm_t norm_i = space->get_norm(i);
            const T mat_ele =  raw_mat_ele * std::conj(charater) * std::sqrt(double(norm_j) / norm_i);
            (columns.contains(state_index) ?  columns[state_index] = mat_ele : columns[state_index] += mat_ele);
        }
    }

    template<typename J, typename T>
    void ref_states(
        const J i,
        const std::unrdered_map<J,T>& row_states, 
        std::unrdered_map<J,T>& rows    
    ) const {
        for(const auto& [state,raw_mat_ele] : row_states){
            const auto& [ref_state,charater] = symmetry.get_refstate(state);
            J j = space->get_index(ref_state);
            const auto norm_j = space->get_norm(j);
            const auto norm_i = space->get_norm(i);
            const T mat_ele =  raw_mat_ele * charater * std::sqrt(double(norm_j) / norm_i);
            (rows.contains(state_index) ? : rows[state_index] = mat_ele : rows[state_index] += mat_ele);
        }
    }

    template<typename J,typename T>
    void calc_rowptr(
        std::vector<operator_string<T>>& pterms,
        std::vector<dense_term<T>>& dterms,
        J rowptr[]
    ) const {

        J n_row = space -> size();

        std::unordered_map<typename bitset_t,T> col_states;
        std::unordered_map<J,T> columns;

        rowptr[0] = 0;
        for(J row = 0;row < n_row;++row)
        {
            col_states.clear();
            columns.clear();
            auto state = space->get_state(row);
            // generate action on states
            for(const auto& pterm : pterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            for(const auto& dterm : dterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            // calculate location of states in basis
            this->ref_state_conj<T>(col_states,columns);
            // insert number of non-zeros elements for this row
            rowptr[row] = columns.size();
        }

        J nnz = 0;
        for(J row = 0;row < n_row;++row){
            J tmp = rowptr[row];
            rowptr[row] = nnz;
            nnz += tmp;
        }
        rowptr[n_row+1] = nnz;

    }

    template<typename J,typename T>
    void calc_matrix(
        std::vector<operator_string<T>>& pterms,
        std::vector<dense_term<T>>& dterms,
        J rowptr[],
        J indices[],
        T values[]
    ) const {

        J n_row = space -> size();

        std::unordered_map<typename bitset_t,T> col_states;
        std::unordered_map<J,T> columns;
        std::vector<std::pair<J,T>> sorted_columns;


        rowptr[0] = 0;
        for(J row = 0;row < n_row;++row)
        {
            col_states.clear();
            columns.clear();
            sorted_columns.clear();

            auto state = space.get_state(row);
            // generate action on states
            for(const auto& pterm : pterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            for(const auto& dterm : dterms){
                pterm.op_dagger<typename bitset_t>(state,col_states);
            }
            // calculate location of states in basis
            this->ref_state_conj<T>(col_states,columns);

            // sort columns
            sorted_columns.insert(columns.begin(), columns.end());
            std::sort(sorted_columns.begin(), sorted_columns.end(), 
                [](std::pair<J,T> lhs,std::pair<J,T> rhs) -> bool 
                {
                    return lhs.first < rhs.first;
                }
            );

            // insert data
            J i = rowptr[row];
            for(const auto& [col,nzval] : sorted_columns){
                indices[i] = col;
                values[i++] = nzval;
            }

        }

    }

    template<typename T, typename X, typename Y>
    void on_the_fly(
        std::vector<operator_string<T>>& pterms,
        std::vector<dense_term<T>>& dterms,
        const Y a,
        const X * x, 
        const Y b, 
        Y  * y
    )const {

        if(b == Y(0.0)){
            std::fill(y,y+space->size(),0)
        }
        else{
            for(size_t index=0;index < space-<size();++index){ y[index] *= b; }
        }

        std::unordered_map<typename bitset_t,T> row_states;
        std::unordered_map<size_t,T> matrix_ele;

        for(J row=0;row < space->size();++row){
            row_states.clear();
            matrix_ele.clear();

            auto state = space.get_state(row);
            // generate action on states
            for(const auto& pterm : pterms){
                pterm.op<typename bitset_t>(state,col_states);
            }
            for(const auto& dterm : dterms){
                pterm.op<typename bitset_t>(state,col_states);
            }
            // calculate location of states in basis
            this->ref_state<T>(row_states,matrix_ele);

            Y total = 0;
            for(const auto& [col,nzval] : matrix_ele){total += nzval * x[col];}
            y[row] += a * total;

            matrix_ele.clear();
        }
    }


};

}

#endif