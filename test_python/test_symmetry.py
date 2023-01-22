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