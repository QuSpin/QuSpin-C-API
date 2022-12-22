#ifndef __QUSPIN_BASIS_TYPES_H__
#define __QUSPIN_BASIS_TYPES_H__

#include "quspin_backend/basis/basis.h"
#include "quspin_backend/basis/symmetry.h"
#include "quspin_backend/matrix/operator.h"

namespace quspin {

template<typename I,typename T>
using dit_symmetry = basis::symmetry<basis::dit_perm<I>, basis::perm_dit<I>, basis::bit_basis::dit_set<I>,T>;

template<typename I,typename J,typename T>
using dit_fullspace_basis = basis::basis<basis::bit_fullspace<I,J>, basis::bit_symmetry<I,T>>

template<typename I,typename J,typename K,typename T>
using dit_subspace_basis = basis::basis<basis::bit_subspace<I,J,K>, basis::bit_symmetry<I,T>>


template<typename I,typename T>
using bit_symmetry = basis::symmetry< basis::bit_perm<I>, basis::perm_bit<I>, basis::bit_basis::bit_set<I>,T>;

template<typename I,typename J,typename T>
using bit_fullspace_basis = basis::basis<basis::bit_fullspace<I,J>, basis::bit_symmetry<I,T>>

template<typename I,typename J,typename K,typename T>
using bit_subspace_basis = basis::basis<basis::bit_subspace<I,J,K>, basis::bit_symmetry<I,T>>

}

#endif