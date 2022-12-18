#ifndef __QUSPIN_BASIS_BASIS_H__
#define __QUSPIN_BASIS_BASIS_H__

template<typename space_t, typename symmetry_t>
class basis
{
private:
    symmetry_t symmetry;
    space_t space;

public:

    typedef space_t::bitset_t bitset_t;
    typedef space_t::index_t index_t;
    typedef space_t::norm_t norm_t;

    basis(symmetry_t _symmetry, space_t _sapce) : symmetry(_symmetry), space(_space) {}
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
        std::unrdered_map<typename bitset_t,T> col_states, 
        std::unrdered_map<typename index_t,T>& columns
    ){
        for(auto const& ele : col_states){
            auto refstate = symmetry.get_refstate(ele.first);
            typename index_t j = state.get_index(refstate.first);
            typename norm_t norm_j = space.get_norm(j);
            typename norm_t norm_i = space.get_norm(i);
            
            columns[state_index] = ele.second * refstate.second * std::sqrt(double(norm_j)/ norm_i);
        }
    }
};


#endif