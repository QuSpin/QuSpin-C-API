#ifndef __QUSPIN_CORE_UTILS_ABI_H__
    #define __QUSPIN_CORE_UTILS_ABI_H__

    #include <quspin/quspin.h>

    namespace quspin_core_abi {
        size_t get_bits(
    const int lhss,
    const int N)
{
    if(lhss < 0 || lhss >= 256){throw std::domain_error("expecting value of lhss to be in range: 0 <= lhss < 256");}
    if(N < 0){throw std::domain_error("expecting value of N to be in range: N>0");}
    const int min_bits = quspin::basis::constants::bits[lhss] * N;
    if(0){}
    else if(min_bits <= 32){return 32;}
    else if(min_bits <= 64){return 64;}
    else if(min_bits <= 128){return 128;}
    else if(min_bits <= 1024){return 1024;}
    else if(min_bits <= 4096){return 4096;}
    else if(min_bits <= 16384){return 16384;}
    else{throw std::runtime_error("The combination of lhss and system size not compatible with QuSpin.");}
}
    }

    #endif