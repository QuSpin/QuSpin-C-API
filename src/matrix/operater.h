#ifndef __QUSPIN_MATRIX_OPERATOR_H__
#define __QUSPIN_MATRIX_OPERATOR_H__

#include <cinttypes>
#include <cmath>
#include <unordered_map>
#include <tuple>
#include <list>
#include <limits>

template<class T>
class dense_term
{
private:
    const int nrow;
    const int nloc;
    const int * loc;
    const T * data;
    const bool * nonzero;

public:
    dense_term(const int _nrow, const T * _data,const bool * _nonzero) : 
    nrow(_nrow), data(_data), nonzero(_nonzero) {}

    ~dense_term(){}

    template<typename I>
    inline void op(I s, std::unordered_map<I,T> &output){
        std::unordered_map<I,T> output();
        int a = get_bits(s,loc,nloc);
        for(int b=0;b<nrow;++b){
            const int i = nrow*a+b;
            if(nonzero[i]){
                I s = Basis::BitBasis::get_sub_bitstring(s,b,loc,nloc);
                output[s] = (output.contains(s) ? output[s] + data[i] : data[i] );
            }
        }
    }

    template<typename I>
    inline void op_dag(I s, std::unordered_map<I,T> &output){
        int a = get_bits(s,loc,nloc);
        for(int b=0;b<nrow;++b){
            const int i = nrow*b+aå;
            if(nonzero[i]){
                I s = Basis::BitBasis::get_sub_bitstring(s,b,loc,nloc);
                output[s] = (output.contains(s) ? output[s] + std::conj(data[i]) : std::conj(data[i]) );
            }
        }
    }

};

template<typename T>
class permutation_term
{
private:
    const int nrow;
    const int nloc;
    const int * loc;
    const int * perm;
    const int * inv_perm;
    const T * data;

public:
    permutation_term(const int _nrow,const int _perm,const int _inv_perm, const T * _data): 
    nrow(_nrow), perm(_perm), inv_perm(_inv_perm), data(_data) {}
    
    ~permutation_term(){}

    template<typename I>
    inline void op(I s, std::unordered_map<I,T> &output){
        int a = Basis::BitBasis::get_sub_bitstring(s,loc,nloc);
        int b = perm[a];
        I s = Basis::BitBasis::set_sub_bitstring(s,b,loc,nloc);
        output[s] = (output.contains(s) ? output[s] + data[i] : data[i] );
    }

    template<typename I>
    inline op_dag(I s, std::unordered_map<I,T> &output){
        int a = Basis::BitBasis::get_sub_bitstring(s,loc,nloc);
        int b = inv_perm[a];
        I s = Basis::BitBasis::set_sub_bitstring(s,b,loc,nloc);
        output[s] = (output.contains(s) ? output[s] + std::conj(data[i]) : std::conj(data[i]) );
    }

};


template<typename T>
class operator
{
private:
    std::vector<permutation_term<T>> pterms;
    std::vector<dense_term<T>> dterms;
    
public:
    operator(/* args */);
    ~operator();

    template<typename Basis,typename J>
    size_t nnz_columns(Basis& basis,J row){
        std::unordered_map<Basis::BitSet,T> output;

        for(auto pterm : pterms){
            pterm->op_dag(basis[row],output);
            
        }

        return 
    }

    template<typename Basis,typename J>
    map<J,T> columns(Basis&,J);
};

operator::operator(/* args */)
{
}

operator::~operator()
{
}



#endif