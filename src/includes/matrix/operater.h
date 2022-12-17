#ifndef __QUSPIN_MATRIX_OPERATOR_H__
#define __QUSPIN_MATRIX_OPERATOR_H__

#include <cinttypes>
#include <cmath>
#include <unordered_map>
#include <tuple>
#include <list>
#include <limits>


namespace quspin {

template<class T>
class dense_term
{
    // project locations onto local hilbert space to evaluate matrix elements
    // good for when the terms have few-body terms where the local matrix fits into cache. 
private: 
    const int lhss; 
    const int nloc;
    const int * loc;
    const T * data;
    const bool * nonzero;

public:
    dense_term(const int _nrow, const T * _data,const bool * _nonzero) : 
    nrow(_nrow), data(_data), nonzero(_nonzero) {}

    ~dense_term(){}

    template<typename I>
    void op(const I s, std::unordered_map<I,T> &output) const {
        const int a = get_bits(s,loc,nloc);
        for(int b=0;b<nrow;++b){
            const int i = nrow*a+b;
            if(nonzero[i]){
                const I r = basis::bit_basis::get_sub_bitstring(s,b,loc,nloc);
                output[r] = (output.contains(r) ? output[r] + data[i] : data[i] );
            }
        }
    }

    template<typename I>
    void op_transpose(const I s, std::unordered_map<I,T> &output) const {
        const int a = get_bits(s,loc,nloc);
        for(int b=0;b<nrow;++b){
            const int i = nrow*b+a;
            if(nonzero[i]){
                const I r = basis::bit_basis::get_sub_bitstring(s,b,loc,nloc);
                output[r] = (output.contains(r) ? output[r] + data[i] : data[i] );
            }
        }
    }

    template<typename I>
    void op_dagger(const I s, std::unordered_map<I,T> &output) const {
        int a = get_bits(s,loc,nloc);
        for(int b=0;b<nrow;++b){
            const int i = nrow*b+a;
            if(nonzero[i]){
                const I r = basis::bit_basis::get_sub_bitstring(s,b,loc,nloc);
                output[r] = (output.contains(r) ? output[r] + std::conj(data[i]) : std::conj(data[i]) );
            }
        }
    }

};

template<typename T>
class operator_string
{
private:
    const int lhss; // local hilbert space size for each term
    const int nloc;  // number of local operators
    const int * loc; // number of local operators in 
    const int * perms; // non-branching operators stored as permutations
    const int * inv_perms; // non-branching operators dagger stored as permutations
    const T * datas; // matrix elements for non-branching operators. 

public:
    operator_string(int _lhss,int _perm,int _inv_perm, T * _data): 
    lhss(_lhss), perm(_perm), inv_perm(_inv_perm), data(_data) {}
    

    ~operator_string(){}

    template<typename I>
    inline void op(const I s, std::unordered_map<I,T> &output) const {
        const int * perm = perms;
        const T * data = datas;
        T m = T(1.0);
        I r = s;
        for(int i=0;i<nloc;++i){
            const int a = basis::bit_basis::get_sub_bitstring(r,loc[i]);
            const int b = perm[a];
            r = basis::bit_basis::set_sub_bitstring(r,a,b,loc[i]);
            m *= data[s_loc];

            if(m==T(0)) break;
            
            // shift to next permutation
            perm += lhss; 
            data += lhss;
        }

        if(m!=T(0)) output[r] = (output.contains(r) ? output[r] + m : m );
    }
    
    template<typename I>
    inline op_transpose(const I s, std::unordered_map<I,T> &output) const {
        const int * perm = inv_perms;
        const T * data = datas;
        T m = T(1.0);
        I r = s;
        for(int i=0;i<nloc;++i){
            const int s_loc = basis::bit_basis::get_sub_bitstring(r,loc[i]);
            const int r_loc = perm[s_loc];
            r = basis::bit_basis::set_sub_bitstring(r,r_loc,loc[i]);
            m *= data[s_loc];

            if(m==T(0)) break;
            
            // shift to next permutation
            perm += lhss; 
            data += lhss;
        }

        if(m!=T(0)) output[r] = (output.contains(r) ? output[r] + m : m );
    }

    template<typename I>
    inline op_dagger(const I s, std::unordered_map<I,T> &output) const {
        const int * perm = inv_perms;
        const T * data = datas;
        T m = T(1.0);
        I r = s;
        for(int i=0;i<nloc;++i){
            const int s_loc = basis::bit_basis::get_sub_bitstring(r,loc[i]);
            const int r_loc = perm[s_loc];
            r = basis::bit_basis::set_sub_bitstring(r,r_loc,loc[i]);
            m *= std::conj(data[s_loc]);

            if(m==T(0)) break;
            
            // shift to next permutation
            perm += lhss; 
            data += lhss;
        }

        if(m!=T(0)) output[r] = (output.contains(r) ? output[r] + m : m );
    }

};



template<typename T>
class operator
{
private:
    std::vector<operator_string<T>> pterms;
    std::vector<dense_term<T>> dterms;

    
public:
    operator(std::vector<operator_string<T>>& pterms,std::vector<dense_term<T>>& dterms){
        self->pterms = pterms;
        self->dterms = dterms;
    }
    ~operator() {}

    template<typename bitset_t>
    void columns(const bitset_t s,unordered_map<bitset_t,T>& output){
        for(auto const& pterm : pterms){
            pterm.op_dagger(s,output);
        }
        for(auto const& dterm : dterms){
            pterm.op_dagger(s,output);
        }
    }

    template<typename bitset_t>
    void rows(const bit_setI s,unordered_map<bitset_t,T>& output){
        for(auto const& pterm : pterms){
            pterm.op(s,output);
        }
        for(auto const& dterm : dterms){
            pterm.op(s,output);
        }
    }

    template<typename basis_t,typename X,typename Y>
    void op(const basis_t &basis,const Y a,const X * x, const V b, Y  * y){
        if(b == Y(0.0)){
            std::fill(y,y+basis.size(),0)
        }


    }





};

}
#endif