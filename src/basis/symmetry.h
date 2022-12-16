#ifndef __QUSPIN_BASIS_SYMMETRY_H__
#define __QUSPIN_BASIS_SYMMETRY_H__


#include <vector>

namespace quspin::basis {

template<typename perm,typename T>
class symmetry
{
private:
    std::vector<perm> perms;
    std::vector<T> chars;

public:
    lattice_symmetry(std::vector<perm> &_perms,std::vector<T> &_chars) : perms(_perms), chars(_chars) {
        assert(_perms.size() == _chars.size());
    }
    ~lattice_symmetry() {}

    size_t size() const {return perms.size();}
    T character(const size_t i) const {return chars[i];}
    perm operator[](const size_t i) const {return perms[i];}

};

}

#endif