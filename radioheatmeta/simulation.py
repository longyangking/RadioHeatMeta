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
        
        
class SimulationPlanar(Simulation):
    pass

class SimulationGrating(Simulation):
    pass

class SimulationPattern(Simulation):
    pass