#ifndef __QUSPIN_BASIS_REFERENCE_H__
#define __QUSPIN_BASIS_REFERENCE_H__

#include <vector>
#include <cmath>
#include <complex>
#include <limits>

#include "utils/functions.h"

namespace quspin::basis {


template<typename dits_or_bits,typename lattice_perm, typename local_perm,typename T>
double check_ref_state(
    const dits_or_bits &s,
    const lattice_symmetry &lat_symm,
    const local_symmetry &loc_symm
){
    double norm=0.0;

    for(int i=0;i<loc_symm.size();++i) 
    {
        const auto r = loc_symm[i].app(s);
        for(int j=0;j<lat_symm.size();++j)
        {
            const auto rr = lat_symm[i].app(r);

            if(rr >  s){return std::numeric_limits<double>::quiet_NaN()};
            if(rr == s) norm +=  std::real(lat_symm.character(j) * loc_symm.character(i));
        }
    }

}

template<typename dits_or_bits,typename lattice_symmetry, typename local_symmetry,typename T>
std::pair<dits_or_bits,T> get_refstate(
    const dits_or_bits s,
    const lattice_symmetry &lat_symm,
    const local_symmetry &loc_symm
){
    
    dits_or_bits ss(s);
    T coeff = T(0);
    
    for(int i=0;i<loc_symm.size();++i) 
    {
        const auto r = loc_symm[i].app(s);
        for(int j=0;j<lat_symm.size();++j)
        {
            const auto rr = lat_symm[i].app(r);

            if(rr > ss){
                ss = rr;
                coeff = lat_symm.character(j) * loc_symm.character(i);
            }
        }
    }

    return make_pair(ss,coeff);

}

}









#endif