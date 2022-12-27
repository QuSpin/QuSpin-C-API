import unittest
from quspin_core.operator import operator_string

class TestBasis(unittest.TestCase):
    def test(self):
        op = operator_string()
        self.assertEqual(op.name,"operator_string")