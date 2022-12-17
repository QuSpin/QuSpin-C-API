#ifndef __QUSPIN_MATRIX_OPERATOR_H__
#define __QUSPIN_MATRIX_OPERATOR_H__

#include <cinttypes>
#include <cmath>
#include <unordered_map>
#include <tuple>
#include <list>
#include <limits>
#include <memory>


namespace quspin {

template<class T>
class dense_term
{
    // project locations onto local hilbert space to evaluate matrix elements
    // good for when the terms have few-body terms where the local matrix fits into cache. 
private: 
    const int lhss; // size of local 
    const int nloc;
    std::unique_ptr<int> loc;
    std::unique_ptr<T> data;
    std::unique_ptr<bool> nonzero;

public:
    dense_term(std::vector<int> _loc,std::vector<std::vector>> _data) : 
    lhss(_data.size()), nloc(_loc.size())
    {
        loc(new int[nloc]);
        data(new T[lhss*lhss]);
        nonzero(new bool[lhss*lhss]);

        // copy to contiguous pointers
        std::copy(_loc.begin(),_loc.end(),loc);

        int ptr = 0; 
        for(int i=0;i<lhss;++i){
            for(int j=0;j<lhss;++j){
                data[ptr] = _data[i][j];
                nonzero[ptr] = (data[ptr] != T(0));
                ++ptr;
            }
        }


    }
    ~dense_term(){}

    template<typename I>
    void op(const I s, std::unordered_map<I,T> &output) const {
        const int a = basis::bit_basis::get_sub_bitstring(s,loc.get(),nloc);
        for(int b=0;b<lhss;++b){ // loop over columns
            const int i = lhss*a+b;
            if(nonzero[i]){
                const I r = basis::bit_basis::set_sub_bitstring(s,b,loc.get(),nloc);
                (output.contains(r) ? output[r] += data[i] : output[r] = data[i] );
            }
        }
    }

    template<typename I>
    void op_transpose(const I s, std::unordered_map<I,T> &output) const {
        const int a = basis::bit_basis::get_sub_bitstring(s,loc.get(),nloc);
        for(int b=0;b<lhss;++b){  // loop over rows
            const int i = lhss*b+a;
            if(nonzero[i]){
                const I r = basis::bit_basis::set_sub_bitstring(s,b,loc.get(),nloc);
                (output.contains(r) ? output[r] += data[i] : output[r] = data[i] );
            }
        }
    }

    template<typename I>
    void op_dagger(const I s, std::unordered_map<I,T> &output) const {
        const int a = basis::bit_basis::get_sub_bitstring(s,loc.get(),nloc);
        for(int b=0;b<lhss;++b){ // loop over rows
            const int i = lhss*b+a;
            if(nonzero[i]){
                const I r = basis::bit_basis::set_sub_bitstring(s,b,loc.get(),nloc);
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
    const int nloc;  // number of local operators
    std::unique_ptr<int> loc; // number of local operators in 
    std::unique_ptr<int> perms; // non-branching operators stored as permutations
    std::unique_ptr<int> inv_perms; // non-branching operators dagger stored as permutations
    std::unique_ptr<T> datas; // matrix elements for non-branching operators. 

public:
    operator_string(std::vector<int> _loc,std::vector<std::vector<int>> _perms, std::vector<std::vector<T>> _datas) : 
    lhss(perms.front().size()), nloc(_locs.size())
    { 

        loc(new int[nloc]);
        perms(new int[nloc*lhss]);
        inv_perms(new int[nloc*lhss]);
        datas(new T[nloc*lhss]);

        int * perm = perms.get();
        int * inv_perm = inv_perms.get();
        T * data = datas.get();

        for(int i=0;i<nloc;i++){
            loc[i] = _loc[i];

            for(int j=0;j<lhss;i++){
                perm[j] = _perms[i][j];
                inv_perm[perm[j]] = j;
                data[j] = _datas[i][j];

            }
            perm += lhss;
            inv_perm += lhss;
            data += lhss;
        }
    }
    

    ~operator_string(){}

    template<typename I>
    inline void op(const I s, std::unordered_map<I,T> &output) const {
        const int * perm = perms.get();
        const T * data = datas.get();
        T m = T(1.0);
        I r = s;
        bool nonzero=true;
        for(int i=0;i<nloc;++i){
            const int a = basis::bit_basis::get_sub_bitstring(r,loc[i]);
            const int b = perm[a];
            r = basis::bit_basis::set_sub_bitstring(r,a,b,loc[i]);
            m *= data[s_loc];

            if(m == T(0)){nonzero=false; break;}
            
            // shift to next permutation
            perm += lhss; 
            data += lhss;
        }

        if( nonzero ) (output.contains(r) ? output[r] += m : output[r] = m );
    }
    
    template<typename I>
    inline op_transpose(const I s, std::unordered_map<I,T> &output) const {
        const int * perm = inv_perms.get();
        const T * data = datas.get();
        T m = T(1.0);
        I r = s;
        bool nonzero = true;
        for(int i=0;i<nloc;++i){
            const int s_loc = basis::bit_basis::get_sub_bitstring(r,loc[i]);
            const int r_loc = perm[s_loc];
            r = basis::bit_basis::set_sub_bitstring(r,r_loc,loc[i]);
            m *= data[s_loc];

            if(m == T(0)){nonzero=false; break;}
            
            // shift to next permutation
            perm += lhss; 
            data += lhss;
        }

        if( nonzero ) (output.contains(r) ? output[r] += m : output[r] = m );
    }

    template<typename I>
    inline op_dagger(const I s, std::unordered_map<I,T> &output) const {
        const int * perm = inv_perms.get();
        const T * data = datas.get();
        T m = T(1.0);
        I r = s;
        bool nonzero=true;
        for(int i=0;i<nloc;++i){
            const int s_loc = basis::bit_basis::get_sub_bitstring(r,loc[i]);
            const int r_loc = perm[s_loc];
            r = basis::bit_basis::set_sub_bitstring(r,r_loc,loc[i]);
            m *= std::conj(data[s_loc]);

            if(m == T(0)){nonzero=false; break;}
            
            // shift to next permutation
            perm += lhss; 
            data += lhss;
        }

        if( nonzero ) (output.contains(r) ? output[r] += m : output[r] = m );
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