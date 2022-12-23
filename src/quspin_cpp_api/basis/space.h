#ifndef __QUSPIN_BASIS_SPACE_H__
#define __QUSPIN_BASIS_SPACE_H__


#include <cmath>
#include <utility>
#include <unordered_map>
#include <algorithm>
// quspin includes
#include "basis/bitbasis/dits.h"
#include "basis/bitbasis/bits.h"


namespace quspin::basis {


// J = {int32, int64}
// I = {uint32, uint64, uint1024, uint16384_t}
// K = {uint8, uint16}


template<typename I,typename J>
class dit_fullspace // sps > 2
{
private:
    const int lhss; // local hilbert space sice
    const J Ns; // number of states
    const I mask; // mask for local hilbert space storage.
    const dit_integer_t bits; // number of bits to store lhss

public:
    typedef dit_set<I> bitset_t;
    typedef J index_t;
    typedef int norm_t;

    dit_fullspace(const int _lhss,const J _Ns) : 
    lhss(_lhss),
    Ns(_Ns), 
    mask(constants::mask[_lhss]), 
    bits(constants::bits[_lhss]) {}

    dit_fullspace(dit_fullspace<I,J>& other) :
    dit_fullspace(other.lhss,other.N,other.Ns,other.mask,other.bits) {}
    
    ~dit_fullspace() {}

    inline J size() const { return Ns;}
    inline J get_Ns() const { return Ns;}
    
    inline bitset_t get_state(const J index) const {
        return ditset<I>(I(Ns-index-1),lhss,mask,bits);
    }

    inline J get_index(const bitset_t& state) const {
        return Ns - integer<J,I>::cast(state.content) - 1;
    }

    static int get_norm(const bitset_t& state) const {return 1;}
    static int get_norm(const J index) const {return 1;}
    
};

template<typename I,typename J,typename K>
class dit_subspace // sps > 2 
{
private:
    const int lhss; // local hilbert space sice
    const I mask; // mask for bits
    const dit_integer_t bits; // number of bits to store lhss

    std::vector<std::pair<I,K>> states;
    std::unordered_map<I,J> index_map;

public:
    typedef dit_set<I> bitset_t;
    typedef J index_t;
    typedef K norm_t;
    
    dit_subspace(const int _lhss, const J Ns_est) : 
    lhss(_lhss), mask(constants::mask[_lhss]), bits(constants::bits[_lhss])
    {
        state.reserve(Ns_est);
        norms.reserve(Ns_est);
        index_map.reserve(Ns_est*2);
    }

    mask(constants::mask[_lhss]), 
    bits(constants::bits[_lhss]) {}
    ~dit_subspace() {}

    inline J size() const { return states.size();}
    inline J get_Ns() const { return states.size();}
    inline size_t nbytes() const {return states.size() * sizeof(std::pair<I,K>)}

    inline bitset_t get_state(const size_t index) const {
        return bitset_t(I(states[index].first),lhss);
    }

    inline J get_index(const bitset_t& state) const {
        return index_map[state.content];
    }

    inline K get_norm(const bitset_t& state) const {
        return states[index_map[state.content]].second;
    }

    inline K get_norm(const J index) const {
        return states[index].second;
    }

    void append(const bitset_t& new_state,const K new_norm){
        if(!index_map.contains(new_state)){
            states.push(make_pair(new_state.content,new_norm))
            index_map[new_state] = states.size();
        }
    }

    void append(const dit_subspace& other){
        for(const auto& [new_state,new_norm] : other.states){
            if(!index_map.contains(new_state)){
                states.push(make_pair(new_state.content,new_norm))
                index_map[new_state] = states.size();
            }
        }
    }

    void sort_states() {
        const bool is_sorted = std::is_sorted(tates.begin(),states.end(),states.begin(),
            [](std::pair<I,L>& lhs,std::pair<I,L>& rhs) -> bool 
            {
                return lhs.first < rhs.first;
            }
        );

        if(!is_sorted){
            std::sort(states.begin(),states.end(),states.begin(),
                [](std::pair<I,L>& lhs,std::pair<I,L>& rhs) -> bool 
                {
                    return lhs.first < rhs.first;
                }
            )

            index_map.clear();
            J index = 0;
            for(const auto& [state,norm] : states){
                index_map[state] = index++;
            }
        }
    }

    void serialize_states(char * output) const {
        char * states_ptr = static_cast<char*>(states.data())
        std::copy(states_ptr,states_ptr+nbytes(),output);
    }

    void deserialize_states(char * input,const size_t nbytes){
        const size_t n_elements = nbytes / sizeof(std::pair<I,K>);
        std::pair<I,K> * input_data = static_cast<std::pair<I,K>*> input;

        states.clear(); // clears current data
        index_map.clear();

        std::copy(input_data, input_data+n_elements,states);

        J index = 0;
        for(const auto& [state,norm] : states){
            index_map[state] = index++;
        }
    }

};


template<typename I,typename J>
class bit_fullspace // sps = 2
{
private:
    const J Ns; // total number of states 
    const dit_integer_t bits; // number of bits to store lhss

public:
    typedef bit_set<I> bitset_t;
    typedef J index_t;
    typedef int norm_t;


    bit_fullspace(const J _Ns) : Ns(_Ns) {}
    ~bit_fullspace() {}

    inline size_t size() const { return Ns;}
    inline size_t get_Ns() const { return Ns;}

   inline bitset_t get_state(const J index) const {
        return ditset<I>(I(Ns-index-1),lhss,mask,bits);
    }

    inline J get_index(const bitset_t& state) const {
        return Ns - integer<J,I>::cast(state.content) - 1;
    }

    static int get_norm(const bitset_t& state) const {return 1;}
    static int get_norm(const J index) const {return 1;}
};

template<typename I,typename J,typename K>
class bit_subspace // sps = 2 
{
private:
    std::vector<std::pair<I,K>> states;
    std::unordered_map<I,J> index_map;
    
public:
    typedef bit_set<I> bitset_t;
    typedef K norm_t;

    bit_subspace(const J Ns_est) {
        state.reserve(Ns_est);
        index_map.reserve(Ns_est*2);
    }
    ~bit_subspace() {}

    inline size_t size() const { return states.size();}
    inline size_t get_Ns() const { return states.size();}
    inline size_t nbytes() const {return states.size() * sizeof(std::pair<I,K>)}

    inline bitset_t get_state(const size_t index) const {
        return bitset_t(I(states[index].first),lhss);
    }

    inline J get_index(const bitset_t& state) const {
        return index_map[state.content];
    }

    inline K get_norm(const bitset_t& state) const {
        return states[index_map[state.content]].second;
    }

    inline K get_norm(const J index) const {
        return states[index].second;
    }

    void append(const bitset_t& new_state,const K new_norm){
        if(!index_map.contains(new_state)){
            states.push(make_pair(new_state.content,new_norm))
            map_index[new_state] = states.size();
        }
    }

    void append(const bit_subspace& other){
        for(const auto& [new_state,new_norm] : other.states){
            if(!index_map.contains(new_state)){
                states.push(make_pair(new_state.content,new_norm))
                index_map[new_state] = states.size();
            }
        }
    }

    void sort_states() {
        const bool is_sorted = std::is_sorted(tates.begin(),states.end(),states.begin(),
            [](std::pair<I,L>& lhs,std::pair<I,L>& rhs) -> bool 
            {
                return lhs.first < rhs.first;
            }
        );

        if(!is_sorted){
            std::sort(states.begin(),states.end(),states.begin(),
                [](std::pair<I,L>& lhs,std::pair<I,L>& rhs) -> bool 
                {
                    return lhs.first < rhs.first;
                }
            )

            index_map.clear();
            J index = 0;
            for(const auto& [state,norm] : states){
                index_map[state] = index++;
            }
        }
    }

    void serialize_states(char * output) const {
        char * states_ptr = static_cast<char*>(states.data())
        std::copy(states_ptr,states_ptr+nbytes(),output);
    }

    void deserialize_states(char * input,const size_t nbytes){
        const size_t n_elements = nbytes / sizeof(std::pair<I,K>);
        std::pair<I,K> * input_data = static_cast<std::pair<I,K>*> input;

        states.clear(); // clears current data
        index_map.clear();

        std::copy(input_data, input_data+n_elements,states);

        J index = 0;
        for(const auto& [state,norm] : states){
            index_map[state] = index++;
        }
    }

};

}
#endif