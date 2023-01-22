from quspin_core._utils import check_is_perm
import unittest

class TestUtils(unittest.TestCase):
    def test_check_is_perm(self):
        self.assertRaises(ValueError,check_is_perm, [1,1,2,3,4])
        self.assertRaises(ValueError,check_is_perm, [0,1,2,3,4,6])
