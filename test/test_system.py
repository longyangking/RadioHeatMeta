import sys
sys.path.append('.')
sys.path.append('..')

import numpy as np 
import unittest
from radioheatmeta.system import Epsilon

class Test_Epsilon(unittest.TestCase):
    def test_init(self):
        try:
            Epsilon(1, "sca")
        except Exception as e:
            self.assertIsNotNone(e)

    def test_epsilon_vals(self):
        epsilon = Epsilon(1, "scalar")
        self.assertTrue(1==epsilon.epsilon_vals)
        epsilon = Epsilon([3,4,5], "diagonal")
        self.assertTrue((np.array([3,4,5])==epsilon.epsilon_vals).all())

class Test_Mur(unittest.TestCase):
    pass

class Test_Lattice(unittest.TestCase):
    pass

if __name__ == "__main__":
    unittest.main()