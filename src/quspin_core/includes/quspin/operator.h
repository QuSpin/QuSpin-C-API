#ifndef __QUSPIN_OPERATOR_H__
#define __QUSPIN_OPERATOR_H__

#include <cinttypes>
#include <cmath>
#include <complex>
#include <unordered_map>
#include <tuple>
#include <list>
#include <limits>
#include <memory>
#include <vector>

#include "quspin/basis/bitbasis/bits.h"
#include "quspin/basis/bitbasis/dits.h"


namespace quspin {

template<class T>
class dense_term
{
    // project locations onto local hilbert space to evaluate matrix elements
    // good for when the terms have few-body terms where the local matrix fits into cache. 
private: 
    const int lhss; // size of local 
    const int nlocs;
    std::vector<int> locs;
    std::vector<T> data;
    std::vector<bool> nonzero;

public:
    dense_term(std::vector<int> _loc,std::vector< std::vector<T> > _data) : 
    lhss(_data.size()), nlocs(_loc.size())
    {

        // copy to contiguous pointers
        std::copy(_loc.begin(),_loc.end(),locs.begin());

        for(auto& d : data){
            for(auto& v : d){
                data.push_back(v);
                nonzero.push_back(v!=T(0));

            }
        }


    }
    ~dense_term(){}

    template<typename bitset_t>
    void op(const bitset_t& s, std::unordered_map<bitset_t,T> &output) const {
        const int a = basis::get_sub_bitstring(s,locs.data(),nlocs);
        for(int b=0;b<lhss;++b){ // loop over columns
            const int i = lhss*a+b;
            if(nonzero[i]){
                const bitset_t r = basis::set_sub_bitstring(s,b,locs.data(),nlocs);
                output[r] = (output.contains(r) ? output[r] + data[i] : data[i] );
            }
        }
    }

    template<typename bitset_t>
    void op_transpose(const bitset_t s, std::unordered_map<bitset_t,T> &output) const {
        const int a = basis::get_sub_bitstring(s,locs.data(),nlocs);
        for(int b=0;b<lhss;++b){  // loop over rows
            const int i = lhss*b+a;
            if(nonzero[i]){
                const bitset_t r = basis::set_sub_bitstring(s,b,locs.data(),nlocs);
                output[r] = (output.contains(r) ? output[r] + data[i] : data[i] );
            }
        }
    }

    template<typename bitset_t>
    void op_dagger(const bitset_t s, std::unordered_map<bitset_t,T> &output) const {
        const int a = basis::get_sub_bitstring(s,locs.data(),nlocs);
        for(int b=0;b<lhss;++b){ // loop over rows
            const int i = lhss*b+a;
            if(nonzero[i]){
                const bitset_t r = basis::set_sub_bitstring(s,b,locs.data(),nlocs);
                (output.contains(r) ? output[r] += std::conj(data[i]) : output[r] = std::conj(data[i]) );
            }
        }
    }

};

template<typename T>
class operator_string
{
private:
    const int lhss; // local hilbert space size for each term
    const int nlocs;  // number of local operators
    std::vector<int> locs; // number of local operators in 
    std::vector<int> perms; // non-branching operators stored as permutations
    std::vector<int> inv_perms; // non-branching operators dagger stored as permutations
    std::vector<T> datas; // matrix elements for non-branching operators. 

public:
    operator_string(std::vector<int> _locs,std::vector<std::vector<int>> _perms, std::vector<std::vector<T>> _datas) : 
    lhss(_perms.front().size()), nlocs(_locs.size())
    { 

        std::copy(_locs.begin(),_locs.end(),locs.begin());

        for(int i=0;i<nlocs;i++){
            datas.insert(datas.end(),_datas[i].begin(),_datas[i].end());
            perms.insert(perms.end(),_perms[i].begin(),_perms[i].end());
            std::vector<int> ip(_perms[i].size());

            for(int j=0;j<_perms[i].size();++j){
                ip[_perms[i][j]] = j;
            }

            inv_perms.insert(inv_perms.begin(),ip.begin(),ip.end());

        }
    }
    
    ~operator_string(){}

    template<typename bitset_t>
    void op(const bitset_t& s, std::unordered_map<bitset_t,T> &output) const {
        const int * perm = perms.data();
        const T * data = datas.data();
        T m = T(1.0);
        bitset_t r(s);
        bool nonzero=true;
        for(int i=0;i<nlocs;++i){
            const int a = quspin::basis::get_sub_bitstring(r,locs[i]);
            const int b = perm[a];
            m *= data[a];
            r = basis::set_sub_bitstring(r,a,b,locs[i]);

            if(m == T(0)){nonzero=false; break;}
            
            // shift to next permutation
            perm += lhss; 
            data += lhss;
        }

        if( nonzero ) output[r] = (output.contains(r) ? output[r] + m : m );
    }
    
    template<typename bitset_t>
    void op_transpose(const bitset_t& s, std::unordered_map<bitset_t,T> &output) const {
        const int * perm = inv_perms.data();
        const T * data = datas.data();
        T m = T(1.0);
        bitset_t r(s);
        bool nonzero = true;
        for(int i=0;i<nlocs;++i){
            const int s_loc = basis::get_sub_bitstring(r,locs[i]);
            const int r_loc = perm[s_loc];
            r = basis::set_sub_bitstring(r,r_loc,locs[i]);
            m *= data[s_loc];

            if(m == T(0)){nonzero=false; break;}
            
            // shift to next permutation
            perm += lhss; 
            data += lhss;
        }

        if( nonzero ) output[r] = (output.contains(r) ? output[r] + m : m );
    }

    template<typename bitset_t>
    void op_dagger(const bitset_t& s, std::unordered_map<bitset_t,T> &output) const {
        const int * perm = inv_perms.data();
        const T * data = datas.data();
        T m = T(1.0);
        bitset_t r(s);
        bool nonzero=true;
        for(int i=0;i<nlocs;++i){
            const int s_loc = basis::get_sub_bitstring(r,locs[i]);
            const int r_loc = perm[s_loc];
            r = basis::set_sub_bitstring(r,r_loc,locs[i]);
            m *= std::conj(data[s_loc]);

            if(m == T(0)){nonzero=false; break;}
            
            // shift to next permutation
            perm += lhss; 
            data += lhss;
        }

        if( nonzero ) output[r] = (output.contains(r) ? output[r] + m : m );
    }

};

}

#ifdef QUSPIN_UNIT_TESTS

namespace quspin{

template class operator_string<double>;
template operator_string<double>::op<bit_set<uint8_t>;
template operator_string<double>::op_dagger<bit_set<uint8_t>;
template operator_string<double>::op_transpose<bit_set<uint8_t>;

template class dense_term<double>;
template dense_term<double>::op<bit_set<uint8_t>;
template dense_term<double>::op_dagger<bit_set<uint8_t>;
template dense_term<double>::op_transpose<bit_set<uint8_t>;

}

#endif

#endif