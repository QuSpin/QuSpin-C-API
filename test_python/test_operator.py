import unittest
from quspin_core.operator import OperatorAPI,OperatorString,NBodyOperator
import numpy as np
from src.quspin_core._utils import _allowed_types

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
        
        locs = np.array([[1,2]])
        perms = np.array([[1,0],[0,1]])
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
        locs = np.array([1,0])
        perms = np.array([[1,0],[0,1]])
        datas = np.array([[1.0,0.5],[0.25,4.254]])

        obj = OperatorString(locs,perms,datas)
        for dtype in _allowed_types.values():
            new_obj = obj.astype(dtype)  
            
            self.assertEqual(new_obj.datas.dtype,dtype)
            self.assertTrue(np.all(datas.astype(dtype)==new_obj.datas))
            self.assertTrue(np.all(locs==new_obj.locs))
            self.assertTrue(np.all(perms==new_obj.perms))
            self.assertEqual(new_obj.lhss,2)
        
    def test_properties(self):
        locs = np.array([1])
        perms = np.array([[1,0]])
        datas = np.array([[1.0,1.0]])

        obj = OperatorString(locs,perms,datas)

        self.assertTrue(np.all(locs==obj.locs))
        self.assertTrue(np.all(perms==obj.perms))
        self.assertTrue(np.all(datas==obj.datas))
        self.assertEqual(obj.lhss,2)
        
class TestNBodyOperator(unittest.TestCase):
    def test_constructor_dtype(self):
        locs = np.array([1,2])
        data = np.arange(4*4).reshape(4,4)
        self.assertRaises(TypeError,NBodyOperator,locs,data)
        
    def test_constructor_data_shape(self):
        locs = np.array([1,2])
        data = np.arange(4*4.0).reshape(4,4,1)
        self.assertRaises(ValueError,NBodyOperator,locs,data)
        
        locs = np.array([1,2])
        data = np.arange(3*4.0).reshape(4,3)
        self.assertRaises(ValueError,NBodyOperator,locs,data)
        
    def test_constructor_locs_shape(self):
        locs = np.array([[1,2]])
        data = np.random.uniform(0,10,size=(4,4))
        self.assertRaises(ValueError,NBodyOperator,locs,data)
        
    def test_constructor_one_body(self):
        locs = np.array([1])
        data = np.random.uniform(0,10,size=(4,4))
        self.assertRaises(ValueError,NBodyOperator,locs,data)
        
    def test_constructor_lhss(self):
        locs = np.array([1,2])
        data = np.random.uniform(0,10,size=(6,6))
        self.assertRaises(ValueError,NBodyOperator,locs,data)
        
    def test_astype(self):
        locs = np.array([1,2])
        data = np.random.uniform(0,10,size=(4,4))
        obj = NBodyOperator(locs,data)

        for dtype in _allowed_types.values():
            new_obj = obj.astype(dtype)  
            
            self.assertEqual(new_obj.data.dtype,dtype)
            self.assertTrue(np.all(data.astype(dtype)==new_obj.data))
            self.assertTrue(np.all(locs==new_obj.locs))
            self.assertEqual(new_obj.lhss,2)
        
        
    def test_properties(self):
        locs = np.array([1,2])
        data = np.random.uniform(0,10,size=(4,4))

        obj = NBodyOperator(locs,data)

        self.assertTrue(np.all(locs==obj.locs))
        self.assertTrue(np.all(data==obj.data))
        self.assertEqual(obj.lhss,2)
        
class TestOperatorAPI(unittest.TestCase):
    def get_valid_strings(self,lhss=2):
        
        p1 = np.random.permutation(np.arange(lhss))
        p2 = np.random.permutation(np.arange(lhss))
        d1 = np.random.normal(size=lhss)
        d2 = np.random.normal(size=lhss)
        op1 = OperatorString([1],p1,d1)
        op2 = OperatorString([2],p2,d2)

        return [op1,op2]
    
    def get_valid_two_body(self):
        sigma_x = np.array([[0.0,1.0],[1.0,0.0]])
        sigma_z = np.array([[-1.0,0.0],[0.0,1.0]])

        op1 = NBodyOperator([0,1],np.kron(sigma_z,sigma_x))
        op2 = NBodyOperator([2,3],np.kron(sigma_z,sigma_x))
        
        return [op1,op2]    
    
    def test_constructor_mo_terms(self):
        self.assertRaises(ValueError,OperatorAPI,[],np.complex128)
    
    def test_constructor_mixed_terms(self):
        self.assertRaises(ValueError,OperatorAPI,[np.eye(4)],np.complex128)
        
        terms = self.get_valid_strings()+self.get_valid_two_body()
        self.assertRaises(ValueError,OperatorAPI,terms,np.complex128)
        
        terms = self.get_valid_strings(2)+self.get_valid_strings(3)
        self.assertRaises(ValueError,OperatorAPI,terms,np.complex128)
        

    def test_constructor(self):
        for dtype in _allowed_types.values():
            a = OperatorAPI(self.get_valid_strings(),dtype)
            b = OperatorAPI(self.get_valid_two_body(),dtype)