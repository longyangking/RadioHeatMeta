import numpy as np 

class Simulation:
    def __init__(self, verbose=False):
        self.__num_G = 0
        self.__omega_list = list()
        self.__num_omega = 0
        self.__target_z = -1

        self.__phi = list()
        self.__kx_start = 0
        self.__kx_end = 0
        self.__ky_start = 0
        self.__ky_end = 0
        self.__num_kx = 0
        self.__num_ky = 0

        self.__lattice = None
        self.__reciprocal_lattice = None
        
        self.__prefactor = 0 
        
        self.__layer_map = None
        self.__material_map = None
        self.__structure = None
    
        self.__target_layer = 0

        self.__Gx_mat = None
        self.__Gy_mat = None
        self.__E_matrices = None
        self.__grand_imaginary_matrices = None
        self.__eps_zz_inv_matrices = None

        self.__source_list = None
        self.__thickness_list = None
        self.__dimension = None
        self.__options = None

    def add_material(self, name, infile):

    def add_material(self, name, omega_list, epsilon_list):

    def set_material(self, name, epsilon, mur, material_type):

    def add_layer(self, name, thickness, material_name):

    def set_layer(self, name, thickness, material_name):

    def set_layer_thickness(self, name, thickness):

    def add_layer_copy(self, name, original_name):

    def delete_layer(self, name):

    def set_source_layer(self, name):

    def set_probe_layer(self, name):

    def set_probe_layer_z_coordinate(self, target_z):

    def set_num_of_G(self, num_G):

    def get_phi(self):

    def get_omega(self):

    def get_epsilon(self, omege_index, position):

    def output_layer_pattern_realization(self, omege_index, name, Nu, Nv, filename):

    def get_num_of_omega(self):

    def init_simulation(self):


        
class SimulationPlanar(Simulation):
    pass

class SimulationGrating(Simulation):
    pass

class SimulationPattern(Simulation):
    pass