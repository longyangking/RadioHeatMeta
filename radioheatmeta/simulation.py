#  Copyright (C) 2020 Yang Long (longyang_123@yeah.net)
# 
#  RadioHeatMeta is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  RadioHeatMeta is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import numpy as np
import system
import rcwa
from fileloader import FileLoader

class Options:
    def __init__(self):
        self.__FMM_rule = "maivfmm"
        self.__integral_method = "gausskronrod"
        self.__polarization = "both"
        self.__print_intermediate = False
        self.__output_flag = ""
        self.__integral_K_parallel = True
        self.__kx_integral_preset = False
        self.__ky_integral_preset = False
        self.__truncation = "circular"

    @property
    def FMM_rule(self):
        return self.__FMM_rule
    
    @FMM_rule.setter
    def FMM_rule(self, rule):
        if rule not in ["maivfmm", "spatialadapative"]:
            raise Exception("unknown FMM rule! only support 2 kinds of rules: maivfmm and spatialadapative")
        self.__FMM_rule = rule

    @property
    def integral_method(self):
        return self.__integral_method

    @integral_method.setter
    def integral_method(self, method):
        if method not in ["gausskronrod", "gausslegendre"]:
            raise Exception("unknown integral method! only support 2 kinds of methods: gausskronrod and gausslegendre")
        self.__integral_method = method

    @property
    def polarization(self):
        return self.__polarization

    @polarization.setter
    def polarization(self, polarization):
        if polarization not in ["te", "tm", "both"]:
            raise Exception("unknown polarization! only support 3 kinds of polarizations: te, tm or both")
        self.__polarization = polarization

    @property
    def output_flag(self):
        return self.__output_flag

    @output_flag.setter
    def output_flag(self, flag):
        self.__output_flag = flag

    @property
    def print_intermediate(self):
        return self.__print_intermediate

    @print_intermediate.setter
    def print_intermediate(self, status):
        self.__print_intermediate = status

    @property
    def integral_K_parallel(self):
        return self.__integral_K_parallel
    
    @integral_K_parallel.setter
    def integral_K_parallel(self, status):
        self.__integral_K_parallel = status

    @property
    def kx_integral_preset(self):
        return self.__kx_integral_preset

    @kx_integral_preset.setter
    def kx_integral_preset(self, status):
        self.__kx_integral_preset = status

    @property
    def ky_integral_preset(self):
        return self.__ky_integral_preset

    @ky_integral_preset.setter
    def ky_integral_preset(self, status):
        self.__ky_integral_preset = status

    @property
    def truncation(self):
        return self.__truncation

    @truncation.setter
    def truncation(self, truncation):
        if truncation not in ["circular", "parallelogramic"]:
            raise Exception("unknown truncation! only support 3 kinds of truncations: circular and parallelogramic")
        self.__truncation = truncation

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
        self.__material_map = None # the list of [name, material class]
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

        self.__num_of_thread = 1
        self.__current_omega_index = -1

        self.verbose = verbose

    def add_material(self, name, infile):
        num_of_material = len(self.__material_map)
        for i in range(num_of_material):
            if self.__material_map[i][0] == name:
                raise Exception("{name}: material already exist!".format(name=name))

        fileloader = FileLoader()
        material = system.Material(name=name,
            omega_list=fileloader.get_omega_list(),
            epsilon_list=fileloader.get_epsilon_list(),
            mur_list=fileloader.get_mur_list()
        )
        if self.verbose:
            print("import material:[{name}] into the simulation".format(name=name))

        self.__material_map.append([name, material])
        self.__structure.add_material(material)


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
            raise Exception("Omega list is None in the simulation!")
        return self.__omega_list

    def get_epsilon(self, omege_index, position):
        MICRON = 1e6 # micro meter unit
        position = MICRON*position # convert the SI unit (m) into the micro meter (um)
        if omege_index < 0:
            raise Exception("Omega index is smaller than zero!")
        if omege_index >= len(self.__omega_list) is None:
            raise Exception("Omega index of of range!")

        # initiate the RCWA for the given omega
        if omege_index != self.__current_omega_index:
            self.__current_omega_index = omege_index
            self.__build_RCWA_matrices()

        layer_index = 0
        offset = 0
        num_of_layer = self.__structure.get_num_of_layer()
        for i in range(num_of_layer):
            if (position[2] > offset) and (position[2] <= offset + self.__thickness_list[i]):
                layer_index = i
                break
            offset += self.__thickness_list[i]

        # only one layer exist
        if (layer_index == 0) and (position[2] > offset):
            layer_index = num_of_layer - 1

        Gx_r, Gx_l = rcwa.mesh_grid(self.__Gx_mat, self.__Gx_mat)
        Gy_r, Gy_l = rcwa.mesh_grid(self.__Gy_mat, self.__Gy_mat)

        Gx_mat = Gx_l - Gx_r
        Gy_mat = Gy_l - Gy_r
        r1, r2, r3, r4 = 0, self.__num_G, self.__num_G, 2*self.__num_G
        pos = (self.__num_G - 1)/2
        if (self.__options.truncation == "circular") and (self.__dimension == "two"):
            pos = 0
        
        # TODO need to think carefully


    def output_layer_pattern_realization(self, omege_index, name, Nu, Nv, filename):
        if (Nu <= 0) or (Nv <= 0):
            raise Exception("Number of point needs to be positive!")

        dx = self.__lattice.bx[0]
        if (Nu > 1):
            dx = dx/(Nu - 1)

        dy = np.hypot(self.__lattice.by[0], self.__lattice.by[1])
        if (Nv > 1):
            dy = dy/(Nv - 1)

        position = np.zeros(3)
        epsilon = np.zeros(9, dtype=complex)

        # TODO Get pattern and output to the file

    def get_num_of_omega(self):
        return len(self.__omega_list)

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

    def __build_RCWA_matrices(self):

    def __reset_simulation(self):

    def __set_target_layer_by_layer(self, layer):

    def __get_structure(self):
    
class SimulationPlanar(Simulation):
    pass

class SimulationGrating(Simulation):
    pass

class SimulationPattern(Simulation):
    pass