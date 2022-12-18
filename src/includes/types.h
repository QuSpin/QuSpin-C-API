#ifndef __QUSPIN_BASIS_TYPES_H__
#define __QUSPIN_BASIS_TYPES_H__

#include "basis/basis.h"
#include "basis/symmetry.h"
#include "matrix/operator.h"

namespace quspin {

template<typename I,typename T>
using dit_symmetry = symmetry<dit_perm<I>,perm_dit<I>,dit_set<I>,T>;

template<typename I,typename J,typename T>
using dit_fullspace_basis = basis<bit_fullspace<I,J>,bit_symmetry<I,T>>

template<typename I,typename J,typename K,typename T>
using dit_subspace_basis = basis<bit_subspace<I,J,K>,bit_symmetry<I,T>>


template<typename I,typename T>
using bit_symmetry = basis::symmetry<bit_perm<I>,perm_bit<I>,bit_set<I>,T>;

template<typename I,typename J,typename T>
using bit_fullspace_basis = basis::basis<bit_fullspace<I,J>,bit_symmetry<I,T>>

template<typename I,typename J,typename K,typename T>
using bit_subspace_basis = basis::basis<bit_subspace<I,J,K>,bit_symmetry<I,T>>

}

#endif