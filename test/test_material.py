import sys
sys.path.append('.')
sys.path.append('..')

import numpy as np 
import unittest
import radioheatmeta.material as Material

class Test_Epsilon(unittest.TestCase):
    def test_init(self):
        try:
            Material.Epsilon(1, "sca")
        except Exception as e:
            self.assertIsNotNone(e)

    def test_epsilon_vals(self):
        epsilon = Material.Epsilon(1, "scalar")
        status = (np.diag(np.array([1,1,1]))==epsilon.epsilon_val).all()
        self.assertTrue(status)

        epsilon = Material.Epsilon([3,4,5], "diagonal")
        status = (np.diag(np.array([3,4,5]))==epsilon.epsilon_val).all()
        self.assertTrue(status)

class Test_Mur(unittest.TestCase):
    def test_init(self):
        try:
            Material.Mu(1, "sca")
        except Exception as e:
            self.assertIsNotNone(e)

    def test_epsilon_vals(self):
        mu = Material.Mu(1, "scalar")
        status = (np.diag(np.array([1,1,1]))==mu.mu_val).all()
        self.assertTrue(status)

        mu = Material.Mu([3,4,5], "diagonal")
        status = (np.diag(np.array([3,4,5]))==mu.mu_val).all()
        self.assertTrue(status)

class Test_Lattice(unittest.TestCase):
    def test_init(self):
        n_omega = 100
        omega_list = np.linspace(1, 10, n_omega)
        epsilon_vals = np.ones(n_omega)
        mu_vals = np.ones(n_omega)
        name = "air"
        material = Material.Material(name=name,
            omega_list=omega_list,
            epsilon_vals=epsilon_vals,
            mu_vals=mu_vals
        )
        self.assertTrue(material.name == name)
        epsilon_type, mu_type = material.get_material_type()
        self.assertTrue(epsilon_type == "scalar")
        self.assertTrue(mu_type == "scalar")

        status = np.zeros(n_omega)
        epsilon_air = Material.Epsilon(1, "scalar")
        mu_air = Material.Mu(1,"scalar")
        for index in range(n_omega):
            status_epsilon = (epsilon_air.epsilon_val == material.get_epsilon_at_index(index)).all()
            status_mu = (mu_air.mu_val == material.get_mu_at_index(index)).all()
            status[index] = status_epsilon and status_mu
        self.assertTrue(status.all())

if __name__ == "__main__":
    unittest.main()