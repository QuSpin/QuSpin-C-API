#ifndef __DIT_PERMS_H__
#define __DIT_PERMS_H__

#include "basis/bitbasis/bit_info.h"
#include "basis/bitbasis/benes_perm.h"

namespace BitBasis {

template <typename I, int D>
class DitPerm
{
private:
    I D_pow[BitInfo<I>::bits];
    const int * _perm;
    const int _length;
public:
    DitPerm(const int * perm,const int _length) : _perm(perm), _length(length) 
    {
        D_pow[0] = 1;
        for(int i=1;i<_length;i++){
            D_pow[i] = I(D) * D_pow[i-1];
        }
    }

    ~DitPerm() {}
};



template <typename I>
class DitPerm<I,2>
{
private:
    Benes::tr_benes<I> benes;
public:
    DitPerm(const int * perm,const int length) 
    {
        Benes::ta_index<I> index;
        for(int i=0;i<BitInfo<I>::bits;i++){index[i] = Benes::no_index;}
        for(int i=0;i<length;i++){index[i] = perm[i];}
        Benes::gen_benes(&benes,index);

    }

    ~DitPerm() {}
};



}


#endif