from quspin_core.symmetry import symmetry_abi
from quspin_core._symmetry_utils import (
    check_bit_lat_args,
    check_dit_lat_args,
    check_bit_loc_args,
    check_dit_loc_args,
    check_bit_args,
    check_dit_args,
    check_args
)
import unittest
import numpy as np



class TestBasis(unittest.TestCase):
    def test_check_bit_lat_args(self):
        # proper arguments
        chars = [1,-1,1]
        args = [
            [0,1,2],
            [1,2,0],
            [2,0,1],
            ]
        
        chars_array = np.array([1,-1,1],dtype=np.complex128)
        
        self.assertEqual(check_bit_lat_args(args,chars),(args,chars_array))
    
    def test_check_dit_lat_args(self):
        pass
    
    def test_check_bit_loc_args(self):
        pass
    
    def test_check_dit_loc_args(self):
        pass
    
    def test_check_bit_args(self):
        pass
    
    def test_check_dit_args(self):
        pass
    
    def test_check_args(self):
        pass