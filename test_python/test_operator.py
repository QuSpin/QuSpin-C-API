import unittest
from quspin_core.operator import OperatorAPI,OperatorString,NBodyOperator
import numpy as np

class TestOperatorString(unittest.TestCase):
    def test_constructor_type(self):
        locs = np.array([1])
        perms = np.array([[1,0]])
        datas = np.array([[1.0,1.0]]).astype(np.int32)
        
        self.assertRaises(TypeError,OperatorString,locs,perms,datas)

    def test_constructor_shape_perms(self):
        locs = np.array([1])
        perms = np.array([[[1,0]]])
        datas = np.array([[1.0,1.0]])
        
        self.assertRaises(ValueError,OperatorString,locs,perms,datas)
        
    def test_constructor_shape_datas(self):
        locs = np.array([1])
        perms = np.array([[1,0]])
        datas = np.array([[[1.0,1.0]]])
        
        self.assertRaises(ValueError,OperatorString,locs,perms,datas)

    def test_constructor_shape_perms_datas(self):
        locs = np.array([[1]])
        perms = np.array([[1,0,1]])
        datas = np.array([[1.0,1.0]])
        
        self.assertRaises(ValueError,OperatorString,locs,perms,datas)
        
    def test_constructor_shape_locs(self):
        locs = np.array([[1]])
        perms = np.array([[1,0]])
        datas = np.array([[1.0,1.0]])
        
        self.assertRaises(ValueError,OperatorString,locs,perms,datas)
        
    def test_constructor_shape_locs_perms(self):
        locs = np.array([1,2])
        perms = np.array([[1,0]])
        datas = np.array([[1.0,1.0]])
        
        self.assertRaises(ValueError,OperatorString,locs,perms,datas)
        
    def test_astype(self):
        locs = np.array([1])
        perms = np.array([[1,0]])
        datas = np.array([[1.0,1.0]])

        obj = OperatorString(locs,perms,datas)
        dtype = np.complex128
        
        obj_complex = obj.astype(dtype)
        self.assertEqual(obj_complex.datas.dtype,dtype)
        
    def test_properties(self):
        locs = np.array([1])
        perms = np.array([[1,0]])
        datas = np.array([[1.0,1.0]])

        obj = OperatorString(locs,perms,datas)

        self.assertTrue(np.all(locs==obj.locs))
        self.assertTrue(np.all(perms==obj.perms))
        self.assertTrue(np.all(datas==obj.datas))
        
class TestNBodyOperator(unittest.TestCase):
    pass