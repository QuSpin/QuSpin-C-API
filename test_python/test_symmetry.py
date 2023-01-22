from quspin_core.symmetry import (
    SymmetryAPI,
    BitPerm,
    PermBit,
    DitPerm,
    PermDit
)
import unittest

class TestBitPerm(unittest.TestCase):

    def test_hash(self):
        perm = [0,1,2,3]
        obj = BitPerm(perm)
        self.assertEqual(hash((2,(0,1,2,3))),hash(obj))
    
    def test_lhss(self):
        perm = [0,1,2,3]
        obj = BitPerm(perm)
        self.assertEqual(obj.lhss,2)
        
class TestPermBit(unittest.TestCase):
    def test_constructor(self):
        obj = PermBit([1,0,1,0])
        self.assertRaises(ValueError,PermBit,[0,1,2,0])

    def test_hash(self):
        obj = PermBit([1,0,1,0])
        self.assertEqual(hash((2,(1,0,1,0))),hash(obj))
    
    def test_lhss(self):
        obj = PermBit([1,0,1,0])
        self.assertEqual(obj.lhss,2)

class TestDitPerm(unittest.TestCase):
    def test_hash(self):
        perm = [0,1,2,3]
        obj = DitPerm(4,perm)
        self.assertEqual(hash((4,(0,1,2,3))),hash(obj))
    
    def test_lhss(self):
        perm = [0,1,2,3]
        obj = DitPerm(4,perm)
        self.assertEqual(obj.lhss,4)
        
class TestPermDit(unittest.TestCase):
    def test_constructor(self):
        args = (3,[0,1],[[0,1,2],[0,1,2],[0,1,2]])
        self.assertRaises(ValueError,PermDit,*args)

        args = (3,[0,1,2],[[0,1,2],[0,1],[0,1,2]])
        self.assertRaises(ValueError,PermDit,*args)
        
    def test_hash(self):
        args = (3,(0,1,2),((0,1,2),(0,1,2),(0,1,2)))
        obj = PermDit(*args)
        self.assertEqual(hash(args),hash(obj))
    
    def test_lhss(self):
        args = (3,(0,1,2),((0,1,2),(0,1,2),(0,1,2)))
        obj = PermDit(*args)
        self.assertEqual(obj.lhss,3)
        
class TestSymmetryAPI(unittest.TestCase):
    def get_functional_args(self):
        L = 4
        bit_lat_args = {}
        dit_lat_args = {}
        for i in range(L):
            perm = [(j+i)%L for j in range(L)]
            
            bit_perm = BitPerm(perm)
            dit_perm = DitPerm(3,perm)
            bit_lat_args[bit_perm] = 1.0
            dit_lat_args[dit_perm] = 1.0
        
        perm_bit_0 = PermBit([0,0,0,0])
        perm_bit_1 = PermBit([0,1,0,1])
        perm_bit_2 = PermBit([1,0,1,0])

        bit_loc_args = {
            perm_bit_0:1.0,
            perm_bit_1:1.0,
            perm_bit_2:1.0,
        }
        
        perm_dit_0 = PermDit(3,[],[])
        perm_dit_1 = PermDit(3,[0,2],[[2,1,0],[2,1,0]])
        perm_dit_2 = PermDit(3,[1,2],[[2,1,0],[2,1,0]])
        
        dit_loc_args = {
            perm_dit_0:1.0,
            perm_dit_1:1.0,
            perm_dit_2:1.0,
        }
        
        return (bit_lat_args,bit_loc_args),(dit_lat_args,dit_loc_args)
        
    
    def test_constructor_type_error(self):
        (bit_lat_args,bit_loc_args), \
        (dit_lat_args,dit_loc_args)  = self.get_functional_args()
        
        symm = next(iter(bit_loc_args.keys()))
        
        bit_lat_args[symm] = 1.0
        args = 2,32,bit_lat_args,bit_loc_args
        self.assertRaises(TypeError,SymmetryAPI,*args)
        
        (bit_lat_args,bit_loc_args), \
        (dit_lat_args,dit_loc_args)  = self.get_functional_args()
        
        symm = next(iter(bit_lat_args.keys()))
        
        bit_loc_args[symm] = 1.0
        args = 2,32,bit_lat_args,bit_loc_args
        self.assertRaises(TypeError,SymmetryAPI,*args)

 
        lat_symm = {BitPerm([0,1,2,3]):1.0}
        loc_symm = {PermDit(2,[],[]):1.0}
        
        args = 2,32,lat_symm,loc_symm
        self.assertRaises(TypeError,SymmetryAPI,*args)    

    
    def test_constructor_value_error(self):
        (bit_lat_args,bit_loc_args), \
        (dit_lat_args,dit_loc_args)  = self.get_functional_args()
        
        symm = DitPerm(4,[0,1,2,3])

        dit_lat_args[symm] = 1.0
        args = 3,32,dit_lat_args,dit_loc_args
        self.assertRaises(ValueError,SymmetryAPI,*args)
        
        (bit_lat_args,bit_loc_args), \
        (dit_lat_args,dit_loc_args)  = self.get_functional_args()   
        
        symm = PermDit(4,[],[])

        dit_loc_args[symm] = 1.0
        args = 3,32,dit_lat_args,dit_loc_args
        self.assertRaises(ValueError,SymmetryAPI,*args)
        
    def test_constructor_no_error(self):
    
        obj = SymmetryAPI(2,32)
        
        (bit_lat_args,bit_loc_args), \
        (dit_lat_args,dit_loc_args)  = self.get_functional_args()

        args = 2,32,bit_lat_args,bit_loc_args
        obj = SymmetryAPI(*args)
        
        args = 3,32,dit_lat_args,dit_loc_args
        obj = SymmetryAPI(*args)
