#ifndef __QUSPIN_BASIS_BASIS_H__
#define __QUSPIN_BASIS_BASIS_H__

#include <unordered_map>
#include <algorithm>
#include <vector>
#include <utility>
#include <memory>

#include "quspin/basis/space.h"
#include "quspin/basis/symmetry.h"
#include "quspin/operator.h"

namespace quspin::basis {

template<typename subspace_t, typename symmetry_t>
class symmetric_basis
{
private:
    symmetry_t symmetry; // OK to copy
    std::shared_ptr<subspace_t> space;

public:

    typedef typename subspace_t::bitset_t bitset_t;
    typedef typename subspace_t::index_t index_t;
    typedef typename subspace_t::norm_t norm_t;

    symmetric_basis(symmetry_t& _symmetry, std::shared_ptr<subspace_t> _space) : symmetry(_symmetry), space(_space) {}
    ~symmetric_basis() {}

    std::shared_ptr<subspace_t> get_space() {return space;}
    symmetry_t get_symmetries() const {return symmetry;}

    template<typename Container,typename Map, typename J>
    void ref_states_conj(
        const J i,
        const Container& col_states, 
              Map& columns
    ) const {
        for(const auto& [state,raw_mat_ele] : col_states){
            const auto& [ref_state,charater] = symmetry.get_refstate(state);
            const J j = space->get_index(ref_state);
            const auto norm_j = space->get_norm(j);
            const auto norm_i = space->get_norm(i);
            const typename Term::value_type mat_ele =  raw_mat_ele * std::conj(charater) * std::sqrt(double(norm_j) / norm_i);
            columns[j] = (columns.contains(j) ?  mat_ele : columns[j] + mat_ele);
        }
    }

    template<typename Container,typename Map, typename J>
    void ref_states(
        const J i,
        const Container& row_states, 
              Map& rows
    ) const {
        for(const auto& [state,raw_mat_ele] : row_states){
            const auto& [ref_state,charater] = symmetry.get_refstate(state);
            J j = space->get_index(ref_state);
            const auto norm_j = space->get_norm(j);
            const auto norm_i = space->get_norm(i);
            const typename Term::value_type mat_ele =  raw_mat_ele * charater * std::sqrt(double(norm_j) / norm_i);
            rows[j] = (rows.contains(j) ? mat_ele : rows[j] + mat_ele);
        }
    }

    template<typename J,typename Term>
    void calc_rowptr(
        std::vector<Term>& terms,
        J rowptr[]
    ) const {

        J n_row = space -> size();


        using bitset_t = typename subspace_t::bitset_t;
        using value_type = typename Term::value_type;

        std::vector<std::pair<bitset_t,value_type>> col_states;
        std::unordered_map<J,value_type> columns;

        col_states.reserve(terms.size());

        rowptr[0] = 0;
        for(J row = 0;row < n_row;++row)
        {
            col_states.clear();
            columns.clear();
            auto state = space->get_state(row);
            // generate action on states
            for(const auto& term : terms){
                term.op_dagger(state,col_states);
            }
            // calculate location of states in basis
            this->ref_states_conj<typename Term::value_type>(col_states,columns);
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

    template<typename J,typename Term>
    void calc_matrix(
        std::vector<Term>& terms,
        Term::value_type values[],
        J rowptr[],
        J indices[]
    ) const {


        using bitset_t = typename subspace_t::bitset_t;
        using value_type = typename Term::value_type;

        std::vector<std::pair<bitset_t,value_type>> col_states;
        std::vector<std::pair<J,value_type>> sorted_columns;
        std::unordered_map<J,value_type> columns;

        col_states.reserve(terms.size());
        sorted_columns.reserve(terms.size());

        J n_row = space -> size();

        rowptr[0] = 0;
        for(J row = 0;row < n_row;++row)
        {
            col_states.clear();
            columns.clear();
            sorted_columns.clear();

            auto state = space.get_state(row);
            // generate action on states
            for(const auto& term : terms){
                pterm.op_dagger(state,col_states);
            }
            // calculate location of states in basis
            this->ref_states_conj<typename Term::value_type>(col_states,columns);

            // sort columns
            sorted_columns.insert(columns.begin(), columns.end());
            std::sort(sorted_columns.begin(), sorted_columns.end(), 
                [](std::pair<J,typename Term::value_type> lhs,std::pair<J,typename Term::value_type> rhs) -> bool 
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

    template<typename Term, typename X, typename Y>
    void on_the_fly(
        std::vector<operator_string<typename Term::value_type>>& terms,
        const Y a,
        const X * x, 
        const Y b, 
        Y  * y
    )const {

        if(b == Y(0.0)){
            std::fill(y,y+space->size(),0);
        }
        else{
            for(typename subspace_t::index_t index=0;index < space->size();++index){ y[index] *= b; }
        }

        using bitset_t = typename subspace_t::bitset_t;
        using value_type = typename Term::value_type;

        std::vector<std::pair<bitset_t,value_type>> col_states;
        std::unordered_map<J,value_type> columns;

        col_states.reserve(terms.size());

        for(typename subspace_t::index_t row=0;row < space->size();++row){
            row_states.clear();
            matrix_ele.clear();

            auto state = space.get_state(row);
            // generate action on states
            for(const auto& term : terms){
                term.op(state,row_states);
            }
            // calculate location of states in basis
            this->ref_states<typename Term::value_type>(row_states,matrix_ele);

            Y total = 0;
            for(const auto& [col,nzval] : matrix_ele){total += nzval * x[col];}
            y[row] += a * total;

            matrix_ele.clear();
        }
    }

};

/*
template<typename subspace_t>
class basis
{
private:
    std::shared_ptr<subspace_t> space;

public:

    typedef typename subspace_t::bitset_t bitset_t;
    typedef typename subspace_t::index_t index_t;
    typedef typename subspace_t::norm_t norm_t;

    basis(std::shared_ptr<subspace_t> _space) : space(_space) {}
    ~basis() {}

    std::shared_ptr<subspace_t> get_space() {return space;}


    template<typename Container, typename J,typename typename Term::value_type>
    void ref_states_conj(
        const J i,
        const Container& col_states, 
        std::unordered_map<J,typename Term::value_type>& columns
    ) const {
        for(const auto& [state,mat_ele] : col_states){
            const J state_index = space->get_index(state);
            columns[state_index] = (columns.contains(state_index) ?  columns[state_index] + std::conj(mat_ele) : std::conj(mat_ele));
        }
    }

    template<typename Container, typename J, typename typename Term::value_type>
    void ref_states(
        const J i,
        const Container& row_states, 
        std::unordered_map<J,typename Term::value_type>& rows    
    ) const {
        for(const auto& [state,mat_ele] : row_states){
            const J state_index = space->get_index(state);
            rows[state_index] = (rows.contains(state_index) ? rows[state_index] + mat_ele : mat_ele);
        }
    }

    template<typename J,typename typename Term::value_type>
    void calc_rowptr(
        std::vector<operator_string<typename Term::value_type>>& pterms,
        std::vector<dense_term<typename Term::value_type>>& dterms,
        J rowptr[]
    ) const {

        J n_row = space -> size();

        std::vector<std::pair<typename subspace_t::bitset_t,typename Term::value_type>> col_states;
        std::unordered_map<J,typename Term::value_type> columns;

        rowptr[0] = 0;
        for(J row = 0;row < n_row;++row)
        {
            col_states.clear();
            columns.clear();
            auto state = space->get_state(row);
            // generate action on states
            for(const auto& pterm : pterms){
                pterm.op_dagger(state,col_states);
            }
            for(const auto& dterm : dterms){
                dterm.op_dagger(state,col_states);
            }
            // calculate location of states in basis
            this->ref_states_conj<typename Term::value_type>(col_states,columns);
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

    template<typename J,typename typename Term::value_type>
    void calc_matrix(
        std::vector<operator_string<typename Term::value_type>>& pterms,
        std::vector<dense_term<typename Term::value_type>>& dterms,
        J rowptr[],
        J indices[],
        typename Term::value_type values[]
    ) const {

        J n_row = space -> size();

        std::vector<std::pair<typename subspace_t::bitset_t,typename Term::value_type>> col_states;
        std::unordered_map<J,typename Term::value_type> columns;
        std::vector<std::pair<J,typename Term::value_type>> sorted_columns;


        rowptr[0] = 0;
        for(J row = 0;row < n_row;++row)
        {
            col_states.clear();
            columns.clear();
            sorted_columns.clear();

            auto state = space.get_state(row);
            // generate action on states
            for(const auto& pterm : pterms){
                pterm.op_dagger(state,col_states);
            }
            for(const auto& dterm : dterms){
                dterm.op_dagger(state,col_states);
            }
            // calculate location of states in basis
            this->ref_states_conj(col_states,columns);

            // sort columns
            sorted_columns.insert(columns.begin(), columns.end());
            std::sort(sorted_columns.begin(), sorted_columns.end(), 
                [](std::pair<J,typename Term::value_type> lhs,std::pair<J,typename Term::value_type> rhs) -> bool 
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

    template<typename typename Term::value_type, typename X, typename Y>
    void on_the_fly(
        std::vector<operator_string<typename Term::value_type>>& pterms,
        std::vector<dense_term<typename Term::value_type>>& dterms,
        const Y a,
        const X * x, 
        const Y b, 
        Y  * y
    )const {

        if(b == Y(0.0)){
            std::fill(y,y+space->size(),0);
        }
        else{
            for(size_t index=0;index < space->size();++index){ y[index] *= b; }
        }

        std::vector<std::pair<typename subspace_t::bitset_t,typename Term::value_type>> row_states;
        std::unordered_map<size_t,typename Term::value_type> matrix_ele;

        for(size_t row=0;row < space->size();++row){
            row_states.clear();
            matrix_ele.clear();

            auto state = space.get_state(row);
            // generate action on states
            for(const auto& pterm : pterms){
                pterm.op(state,row_states);
            }
            for(const auto& dterm : dterms){
                dterm.op(state,row_states);
            }
            // calculate location of states in basis
            this->ref_state(row_states,matrix_ele);

            Y total = 0;
            for(const auto& [col,nzval] : matrix_ele){total += nzval * x[col];}
            y[row] += a * total;

            matrix_ele.clear();
        }
    }


};

*/

}// end namespace

#ifdef QUSPIN_UNIT_TESTS


// TODO: instaiate code here

#endif

#endif