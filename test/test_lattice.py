import sys
sys.path.append('.')
sys.path.append('..')

import numpy as np 
import unittest
from radioheatmeta.lattice import Lattice

class Test_Lattice(unittest.TestCase):
    def test_init(self):
        try:
            Lattice("a")
        except Exception as e:
            self.assertIsNotNone(e)

        try:
            Lattice("two",ax="sa")
        except Exception as e:
            self.assertIsNotNone(e)

    def test_Lattice_G(self):
        lattice = Lattice("two", ax=[1,0])
        self.assertTrue(lattice.get_dimension() == "two")

        lattice.init_lattice(100)
        nG = lattice.get_nG()
        self.assertTrue((nG <= 100) and (nG >= 1))

        lattice.init_lattice(100, truncation="circular")
        nG = lattice.get_nG()
        self.assertTrue((nG <= 100) and (nG >= 1))

if __name__ == "__main__":
    unittest.main()
    
    #lattice = Lattice("two", ax=[1,0], verbose=True)
    #lattice.init_lattice(100, truncation="circular")