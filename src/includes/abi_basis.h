#ifndef __QUSPIN_ABI_BASIS_H__
#define __QUSPIN_ABI_BASIS_H__

// basis type templates I,J,K,T
// basis dispatching only in I:
// I = uint32_t | uint64_t | uint1024_t | uint4096_t | uint16384_t
// J = npy_intp | npy_intp | npy_intp   | npy_intp   | npy_intp
// K = uint8_t  | uint8_t  | int        | int        | int
// T = c128     | c128     | c128       | c128       | c128

// basis type templates J,T
// matrix dispatching on all combinations:
// J = {int32_t, int64_t}
// T = {int8_t, float32_t, float64_t, complex64_t, complex128_t}



class symmetric_bit_basis
{


private:
    unqie_ptr<void> data;
    const bool symmetric;

public:
    NPY_TYPES I_dtype

    bit_basis() {

    }
    ~bit_basis(){}    

};

#endif