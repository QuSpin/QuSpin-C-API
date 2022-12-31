#ifndef __QUSPIN_OPERATOR_H__
#define __QUSPIN_OPERATOR_H__

#include <cinttypes>
#include <cmath>
#include <complex>
#include <array>
#include <unordered_map>
#include <utility>
#include <list>
#include <limits>
#include <memory>
#include <vector>

#include "quspin/basis/bitbasis/bits.h"
#include "quspin/basis/bitbasis/dits.h"
#include "quspin/utils/functions.h"


namespace quspin {




template<class T,std::size_t N>
class N_body_dits
{
    //
private: 
    basis::dit_integer_t lhss; // size of local 
    int dim;
    std::array<int,N> locs;
    std::vector<T> data;
    std::vector<bool> nonzero;

public:
    typedef T value_type;

    N_body_dits(const basis::dit_integer_t _lhss,std::vector<int> _locs,std::vector<T> &_data) : 
    lhss(_lhss)
    {
        dim = 1; for(int i=0;i<N;i++){dim *= _lhss;}

        assert(_data.size() == dim*dim);
        assert(_locs.size() == N);
        // copy to contiguous pointers
        data.insert(_data.begin(),_data.end(),data.begin());
        std::copy(_locs.begin(),_locs.end(),locs.begin());
        
        std::transform(data.begin(),data.end(),nonzero.begin(),
            [](const T& value) -> bool 
            {
                return value != T(0);
            }
        );

    }
    ~N_body_dits(){}

    template<typename bitset_t>
    void op(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        const int a = basis::get_sub_bitstring(s,locs);
        for(int b=0;b<dim;++b){ // loop over columns
            const int i = dim*a+b;
            if(nonzero[i]){
                const bitset_t r = basis::set_sub_bitstring(s,b,locs);
                output.push_back(std::make_pair(r,data[i]));
            }
        }
    }

    template<typename bitset_t>
    void op_transpose(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        const int a = basis::get_sub_bitstring(s,locs);
        for(int b=0;b<dim;++b){  // loop over rows
            const int i = dim*b+a;
            if(nonzero[i]){
                const bitset_t r = basis::set_sub_bitstring(s,b,locs);
                output.push_back(std::make_pair(r,data[i]));
            }
        }
    }

    template<typename bitset_t>
    void op_dagger(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        const int a = basis::get_sub_bitstring(s,locs);
        for(int b=0;b<dim;++b){ // loop over rows
            const int i = dim*b+a;
            if(nonzero[i]){
                const bitset_t r = basis::set_sub_bitstring(s,b,locs);
                output.push_back(std::make_pair(r,conj(data[i])));
            }
        }
    }

};

template<class T,std::size_t N>
class N_body_bits
{
    //
private: 
    enum {dim = integer_pow<N,2>::value};

    std::array<int,N> locs;
    std::array<T,dim*dim> data;
    std::array<bool,dim*dim> nonzero;

public:
    typedef T value_type;

    N_body_bits(std::vector<int> _locs,std::vector<T> &_data)
    {
        assert(_data.size() == dim*dim);
        assert(_locs.size() == N);
        // copy to contiguous pointers
        std::copy(_data.begin(),_data.end(),data.begin());
        std::copy(_locs.begin(),_locs.end(),locs.begin());
        
        std::transform(data.begin(),data.end(),nonzero.begin(),
            [](const T& value) -> bool 
            {
                return value != T(0);
            }
        );

    }
    ~N_body_bits(){}

    template<typename bitset_t>
    void op(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        const int a = basis::get_sub_bitstring(s,locs);

        for(int b=0;b<dim;++b){ // loop over columns
            const int i = dim*a+b;
            if(nonzero[i]){
                bitset_t r = basis::set_sub_bitstring(s,b,locs);
                output.push_back(std::make_pair(r,data[i]));
            }
        }
    }

    template<typename bitset_t>
    void op_transpose(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        const int a = basis::get_sub_bitstring(s,locs);

        for(int b=0;b<dim;++b){  // loop over rows
            const int i = dim*b+a;
            if(nonzero[i]){
                bitset_t r = basis::set_sub_bitstring(s,b,locs);
                output.push_back(std::make_pair(r,data[i]));
            }
        }
    }

    template<typename bitset_t>
    void op_dagger(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        const int a = basis::get_sub_bitstring(s,locs);

        for(int b=0;b<dim;++b){ // loop over rows
            const int i = dim*b+a;
            if(nonzero[i]){
                bitset_t r = basis::set_sub_bitstring(s,b,locs);
                output.push_back(std::make_pair(r,conj(data[i])));
            }
        }
    }

};

template<typename T>
class operator_string // generic operator
{
private:
    const int lhss; // local hilbert space size for each term
    const int nlocs;  // number of local operators
    std::vector<int> locs; // number of local operators in 
    std::vector<int> perms; // non-branching operators stored as permutations
    std::vector<int> inv_perms; // non-branching operators dagger stored as permutations
    std::vector<T> datas; // matrix elements for non-branching operators. 

public:

    typedef T value_type;

    operator_string(std::vector<int> _locs,std::vector<std::vector<int>> _perms, std::vector<std::vector<T>> _datas) : 
    lhss(_perms.front().size()), nlocs(_locs.size())
    { 

        locs.insert(locs.end(),_locs.begin(),_locs.end());

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
    void op(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        const int * perm = perms.data();
        const T * data = datas.data();
        T m = T(1.0);
        bitset_t r(s);
        bool nonzero=true;
        for(int i=0;i<nlocs;++i){
            const int a = quspin::basis::get_sub_bitstring(r,locs[i]);
            const int b = perm[a];
            m *= data[a];
            r = basis::set_sub_bitstring(r,b,locs[i]);

            if(m == T(0)){nonzero=false; break;}
            
            // shift to next permutation
            perm += lhss; 
            data += lhss;
        }

        if( nonzero ) output.push_back(std::make_pair(r,m));
    }
    
    template<typename bitset_t>
    void op_transpose(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
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

        if( nonzero ) output.push_back(std::make_pair(r,m));
    }

    template<typename bitset_t>
    void op_dagger(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        const int * perm = inv_perms.data();
        const T * data = datas.data();
        T m = T(1.0);
        bitset_t r(s);
        bool nonzero=true;
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

        if( nonzero ) output.push_back(std::make_pair(r,conj(m)));
    }

};

}

#ifdef QUSPIN_UNIT_TESTS

namespace quspin {

typedef basis::bit_set<uint8_t> bs;

template class operator_string<double>;
template void operator_string<double>::op<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;
template void operator_string<double>::op_dagger<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;
template void operator_string<double>::op_transpose<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;

template class N_body_dits<double,2>;
template void N_body_dits<double,2>::op<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;
template void N_body_dits<double,2>::op_dagger<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;
template void N_body_dits<double,2>::op_transpose<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;


template class N_body_bits<double,2>;
template void N_body_bits<double,2>::op<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;
template void N_body_bits<double,2>::op_dagger<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;
template void N_body_bits<double,2>::op_transpose<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;



}

TEST_CASE("operator_string.op"){
    using namespace quspin;



    // SzSz
    std::vector<std::vector<double>> datas = {{-0.5,0.5},{-0.5,0.5}};
    std::vector<std::vector<int>> perms = {{0,1},{0,1}};
    std::vector<int> locs = {0,1};

    operator_string<double> SzSz(locs,perms,datas);

    // // SxSx
    // datas = {{0.5,0.5},{0.5,0.5}};
    // perms = {{1,0},{1,0}};
    // locs = {0,1};

    // operator_string<double> SxSx(locs,datas,perms);


    // // SySy
    // datas = {{-0.5,0.5},{0.5,-0.5}};
    // perms = {{1,0},{1,0}};
    // locs = {0,1};

    // operator_string<double> SySy(locs,datas,perms);


    
}

TEST_CASE("operator_string.op_dagger"){
    
}

TEST_CASE("operator_string.op_transpose"){
    
}

TEST_CASE("N_body_bits<double,2>.op"){

}

TEST_CASE("N_body_bits<double,2>.op_dagger"){
    
}

TEST_CASE("N_body_bits<double,2>.op_transpose"){
    
}

TEST_CASE("N_body_dits<double,2>.op"){

}

TEST_CASE("N_body_dits<double,2>.op_dagger"){
    
}

TEST_CASE("N_body_dits<double,2>.op_transpose"){
    
}

#endif

#endif