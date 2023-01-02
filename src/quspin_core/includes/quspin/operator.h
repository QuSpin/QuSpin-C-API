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
    std::vector<T> inv_datas; // matrix elements for dagger operator.

public:

    typedef T value_type;

    operator_string(std::vector<int> _locs,std::vector<std::vector<int>> _perms, std::vector<std::vector<T>> _datas) : 
    lhss(_perms.front().size()), nlocs(_locs.size())
    { 

        assert(_locs.size() == _perms.size());
        assert(_locs.size() == _datas.size());

        locs.insert(locs.end(),_locs.begin(),_locs.end());

        for(int i=0;i<_locs.size();i++){

            std::vector<int> perm(_perms[i].begin(),_perms[i].end()),inv_perm(_perms[i].size());
            std::vector<T> data(_datas[i].begin(),_datas[i].end()),inv_data(_datas[i].size());

            assert(perm.size() == lhss);
            assert(data.size() == lhss);

            int j = 0;
            for(const int p : perm){
                assert((p >= 0) && (p < lhss));
                inv_perm[p] = j;
                inv_data[p] = data[j];
                j++;
            }

            perms.insert(perms.end(),perm.begin(),perm.end());
            inv_perms.insert(inv_perms.end(),inv_perm.begin(),inv_perm.end());

            datas.insert(datas.end(),data.begin(),data.end());
            inv_datas.insert(inv_datas.end(),inv_data.begin(),inv_data.end());

        }
    }
    
    ~operator_string(){}

    template<typename bitset_t>
    void op(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        const int * perm = perms.data();
        const T * data = datas.data();
        T m = T(1.0);
        bitset_t r(s);

        bool nonzero = true;
        size_t ptr =0;
        
        for(int i=0;i<nlocs;++i){
            const int s_loc = quspin::basis::get_sub_bitstring(r,locs[i]);
            const size_t ind = ptr + s_loc;

            m *= datas[ind];
            r = basis::set_sub_bitstring(r,perms[ind],locs[i]);

            if(m == T(0)){nonzero=false; break;}
            
            // shift to next permutation
            ptr += lhss;
        }

        if( nonzero ) output.push_back(std::make_pair(r,m));
    }
    
    template<typename bitset_t>
    void op_transpose(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        T m = T(1.0);
        bitset_t r(s);
        
        bool nonzero = true;
        size_t ptr = inv_perms.size()-lhss;

        for(int i=nlocs-1;i > -1;--i){
            const int s_loc = basis::get_sub_bitstring(r,locs[i]);
            const size_t ind = ptr + s_loc;

            m *= inv_datas[ind];
            r = basis::set_sub_bitstring(r,inv_perms[ind],locs[i]);

            if(m == T(0)){nonzero=false; break;}
            ptr -= lhss;
        }

        if( nonzero ) output.push_back(std::make_pair(r,m));
    }

    template<typename bitset_t>
    void op_dagger(const bitset_t& s, std::vector<std::pair<bitset_t,T>> &output) const {
        T m = T(1.0);
        bitset_t r(s);
        
        bool nonzero = true;
        size_t ptr = inv_perms.size()-lhss;

        for(int i=nlocs-1;i > -1;--i){
            const int s_loc = basis::get_sub_bitstring(r,locs[i]);
            const size_t ind = ptr + s_loc;

            m *= inv_datas[ind];
            r = basis::set_sub_bitstring(r,inv_perms[ind],locs[i]);

            if(m == T(0)){nonzero=false; break;}
            
            ptr -= lhss;
        }

        if( nonzero ) output.push_back(std::make_pair(r,conj(m)));
    }

};

}

#ifdef QUSPIN_UNIT_TESTS

namespace quspin { // test cases

typedef basis::bit_set<uint8_t> bs;
typedef basis::dit_set<uint8_t> ds;

// qubits

template class operator_string<double>;
template void operator_string<double>::op<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;
template void operator_string<double>::op_dagger<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;
template void operator_string<double>::op_transpose<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;

template class N_body_bits<double,2>;
template void N_body_bits<double,2>::op<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;
template void N_body_bits<double,2>::op_dagger<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;
template void N_body_bits<double,2>::op_transpose<bs>(const bs& , std::vector<std::pair<bs,double>>&) const;

// qudits

template class N_body_dits<double,2>;
template void N_body_dits<double,2>::op<ds>(const ds& , std::vector<std::pair<ds,double>>&) const;
template void N_body_dits<double,2>::op_dagger<ds>(const ds& , std::vector<std::pair<ds,double>>&) const;
template void N_body_dits<double,2>::op_transpose<ds>(const ds& , std::vector<std::pair<ds,double>>&) const;

}



TEST_CASE("operator_string.op"){
    using namespace quspin;
    operator_string<double> * H;

    // qudits
    std::vector<std::pair<ds,double>> dit_output;
    basis::dit_set<uint8_t> dit_state({0,1,2,1},3);

    H = new operator_string<double>({0,1},{{0,1,2},{0,1,2}}, {{-1.0,0.0,1.0},{-1.0,0.0,1.0}});

    H->op(dit_state,dit_output);
    CHECK(dit_output.size() == 0);

    delete H;

    H = new operator_string<double>({0,2},{{0,1,2},{0,1,2}}, {{-1.0,0.0,1.0},{-1.0,0.0,1.0}});

    dit_output.clear();
    H->op(dit_state,dit_output);
    CHECK(dit_output.size() == 1);
    CHECK(dit_output[0].first.to_string() ==  "0 1 2 1 ");
    CHECK(dit_output[0].second == -1.0);

    delete H;

    H = new operator_string<double>({0,2,3},{{0,1,2},{1,2,0},{2,0,1}}, {{1.0,2.0,3.0},{1.0,2.0,3.0},{1.0,2.0,3.0}});

    dit_output.clear();
    H->op(dit_state,dit_output);
    CHECK(dit_output.size() == 1);
    CHECK(dit_output[0].first.to_string() ==  "0 1 0 0 ");
    CHECK(dit_output[0].second == 6.0);

    delete H;


    // qubits

    std::vector<std::pair<bs,double>> bit_output;
    basis::bit_set<uint8_t> bit_state({0,1,1,0,1,0,0,1});

    // \sigma_z \sigma_z
    H = new operator_string<double>({0,1},{{0,1},{0,1}}, {{-1.0,1.0},{-1.0,1.0}});
    H->op(bit_state,bit_output);

    CHECK(bit_output.size() == 1);
    CHECK(bit_output[0].first.to_string() == "0 1 1 0 1 0 0 1 ");
    CHECK(bit_output[0].second == -1.0);

    delete H;

    H = new operator_string<double>({0,3},{{0,1},{0,1}}, {{-1.0,1.0},{-1.0,1.0}});
    H->op(bit_state,bit_output);

    CHECK(bit_output.size() == 2);
    CHECK(bit_output[1].first.to_string() == "0 1 1 0 1 0 0 1 ");
    CHECK(bit_output[1].second == 1.0);


    H = new operator_string<double>({0,1,4},{{0,1},{0,1},{1,0}}, {{-1.0,1.0},{-1.0,1.0},{0.5,1.0}});
    bit_output.clear();
    H->op(bit_state,bit_output);

    CHECK(bit_output.size() == 1);
    CHECK(bit_output[0].first.to_string() ==  "0 1 1 0 0 0 0 1 ");
    CHECK(bit_output[0].second == -1.0);

    delete H;

    H = new operator_string<double>({0,1,0},{{0,1},{0,1},{1,0}}, {{-1.0,1.0},{-1.0,1.0},{0.5,1.0}});
    bit_output.clear();
    H->op(bit_state,bit_output);

    CHECK(bit_output.size() == 1);
    CHECK(bit_output[0].first.to_string() ==  "1 1 1 0 1 0 0 1 ");
    CHECK(bit_output[0].second == -0.5);

    delete H;


}


TEST_CASE("operator_string.op_dagger"){
 
}

TEST_CASE("operator_string.op_transpose"){
    using namespace quspin;

    typedef double T;
    operator_string<T> * H;
    std::vector<std::pair<ds,T>> output;
    basis::dit_set<uint8_t> state({2,10},16);

    std::vector<int> perm;
    std::vector<T> data;

    for(int i=0;i<16;i++){
        int ip = (i+1)%16;
        perm.push_back(ip);
        data.push_back(std::sqrt(ip));
    }

    H = new operator_string<T>({0,1},{perm,perm},{data,data});

    H->op_transpose(state,output);
    CHECK(output.size() == 1);
    CHECK(output[0].first.to_string() ==  "1 9 ");
    REQUIRE(output[0].second == doctest::Approx(std::sqrt(2*10)));

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