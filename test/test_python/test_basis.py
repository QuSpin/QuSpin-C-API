import unittest
from quspin_core.basis import bit_basis

class TestBasis(unittest.TestCase):
    def test(self):
        b = bit_basis()
        self.assertEqual(b.name=="operator")