import numpy as np 

class Simulation:
    def __init__(self, verbose=False):
        self.__num_G = 0
        self.__omega_list = None
        self.__num_omega = 0
        self.__target_z = -1

        self.__phi = None
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
        self.__dimension = "no"
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
        if self.__phi is None:
            raise Exception("Phi is None, Please integrate firstly!")
        return self.__phi

    def get_omega(self):
        if self.__omega_list is None:
            raise Exception("Omega list is None!")
        return self.__omega_list

    def get_epsilon(self, omege_index, position):

    def output_layer_pattern_realization(self, omege_index, name, Nu, Nv, filename):

    def get_num_of_omega(self):

    def init_simulation(self):

    def get_phi_at_kx_ky(self, omega_index, kx, ky=0):

    def get_num_of_G(self):

    def output_sys_info(self):

    def opt_print_intermediate(self, output_flat=""):

    def opt_only_compute_TE(self):

    def opt_only_compute_TM(self):

    def opt_set_lattice_truncation(self, truncation):

    def set_thread(self, num_thread):

    def set_kx_integral(self, points, end=0):

    def set_kx_integral_sym(self, points, end=0):

    def set_ky_integral(self, points, end=0):

    def set_ky_integral_sym(self, points, end=0):

    def integrate_ky_ky(self):

    def integrate_kx_ky_MPI(self, rank, size):

    def integrate_kx_ky_GPU(self):

    def __integrate_kx_ky_internal(self, start, end, parallel=False, rank=0):

    def build_RCWA_matrices(self):

    def reset_simulation(self):

    def set_target_layer_by_layer(self, layer):

    def get_structure(self):

        
class SimulationPlanar(Simulation):
    pass

class SimulationGrating(Simulation):
    pass

class SimulationPattern(Simulation):
    pass