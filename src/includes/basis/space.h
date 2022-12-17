#ifndef __QUSPIN_BASIS_SPACE_H__
#define __QUSPIN_BASIS_SPACE_H__


#include <cmath>
#include "basis/bitbasis/dits.h"
#include "basis/bitbasis/bits.h"


namespace quspin::basis {


// J = {int32, int64}
// I = {uint32, uint64, uint1024, uint16384_t}
// K = {uint8, uint16}

template<typename J>
J calc_Ns_fullspace(const int lhss,const int N)
{
    J Ns = 1;
    for(int i=0;i<N;i++){Ns *= lhss;}
    return Ns;
}

template<typename I,typename J>
class dit_fullspace // sps > 2
{
private:
    const int lhss; // local hilbert space sice
    const int N; // number of lattice sites
    const J Ns; // number of states
    const I mask; // mask for local hilbert space storage.
    const dit_integer_t nbits; // number of bits to store lhss

public:
    typedef I bitset_t;
    typedef J index_t;

    dit_fullspace(const int _lhss,const int _N) : 
    lhss(_lhss), N(_N) 
    Ns(calc_Ns_fullspace<J>(_lhss,_N)), 
    mask(constants::mask[_lhss]), 
    nbits(constants::nbits[_lhss]) {}

    ~dit_fullspace() {}

    inline J size() const { return Ns;}
    inline J get_Ns() const { return Ns;}
    inline int get_N() const {return N;}
    inline dit_set<I> operator[](const J index) const {
        return ditset<I>(I(Ns-index-1),lhss,mask,nbits);
    }
    inline J operator[](const dit_set<I> state) const {
        return Ns - integer_cast<J,I>(state.content) - 1;
    }
};


template<typename I,typename J,typename K>
class dit_subspace // sps > 2 
{
private:
    const int lhss; // local hilbert space sice
    const int N; // number of lattice sites
    const I mask; // mask for nbits
    const dit_integer_t nbits; // number of bits to store lhss
    std::vector<I> states;
    std::vector<K> norms;
    std::unordered_map<I,J> index_map;


public:
    typedef I bitset_t;
    typedef J index_t;

    dit_fullspace(const int _lhss,const int N) : 
    lhss(_lhss), N(_N),
    mask(constants::mask[_lhss]), 
    nbits(constants::nbits[_lhss]) {}
    ~dit_fullspace() {}

    inline J size() const { return states.size();}
    inline J get_Ns() const { return states.size();}
    inline int get_N() const { return N;}

    dit_set<I> operator[](const size_t index) const {
        return dit_set<I>(I(states[index]),lhss);
    }

    J operator[](const dit_set<I> state) const {
        return index_map[state.content];
    }


};


template<typename I,typename J>
class bit_fullspace // sps = 2
{
private:
    const int N; // number of lattice sites
    const J Ns; // total number of states 
    const dit_integer_t nbits; // number of bits to store lhss

public:
    typedef I bitset_t;
    typedef J index_t;


    bit_fullspace(const int _N) : 
    N(_N), Ns(calc_Ns_fullspace<J>(2,_N)) {}

    ~bit_fullspace() {}

    inline size_t size() const { return Ns;}
    inline size_t get_Ns() const { return Ns;}
    inline int get_N() const { return N;}

    inline bit_set<I> operator[](const J index) const {
        return bit_set<I>(I(Ns-index-1));
    }
    inline J operator[](const bit_set<I> state) const {
        return Ns - integer_cast<J,I>(state.content) - 1;
    }
};

template<typename I,typename J,typename K>
class bit_subspace // sps = 2 
{
private:
    const int N; // number of lattice sites
    std::vector<I> states;
    std::vector<K> norms;
    std::unordered_map<I,J> index_map;


public:
    typedef I bitset_t;
    typedef J index_t;

    bit_subspace(const int N) : 
    N(_N),
    ~bit_subspace() {}

    inline size_t size() const { return states.size();}
    inline size_t get_Ns() const { return states.size();}
    inline int get_N() const { return N;}

    bit_set<I> operator[](const size_t index) const {
        return bit_set<I>(I(states[index]),lhss);
    }

    J operator[](const bit_set<I> state) const {
        return index_map[state.content];
    }

    template<typename gen_t,typename lattice_symmetry, typename local_symmetry> 
    void grow_basis(gen_t generator,lattice_symmetry &lat_symm, local_symmetry &loc_symm){
        for(const auto& s : generator){
            
        }
    }
};

}
#endif