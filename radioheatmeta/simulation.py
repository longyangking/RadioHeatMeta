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
from radioheatmeta.calculation import RCWA

class Options:
    def __init__(self):
        self.__polarization = "both"
        self.verbose = False
        self.__truncation = "circular"
        self.nG = 100

    @property
    def polarization(self):
        return self.__polarization

    @polarization.setter
    def polarization(self, polarization):
        if polarization not in ["te", "tm", "both"]:
            raise Exception("unknown polarization! only support 3 kinds of polarizations: te, tm or both")
        self.__polarization = polarization

    @property
    def truncation(self):
        return self.__truncation

    @truncation.setter
    def truncation(self, truncation):
        if truncation not in ["circular", "parallelogramic"]:
            raise Exception("unknown truncation! only support 2 kinds of truncations: circular and parallelogramic")
        self.__truncation = truncation

class Simulation:
    def __init__(self, structure=None, options=None, verbose=False):
        self.__structure = structure
        self.__omega_list = structure.get_omega_list()

        self.__options = options
        if options is None:
            self.__options = Options()
        self.verbose = verbose

        self.__target_z = 0
        self.__target_layer_index = 0
    
    def init_simulation(self):
        pass

    def add_layer(self, layer):
        self.__structure.add_layer(layer)

    def set_layer(self, index, layer):
        self.__structure.set_layer_by_index(
            index=index,
            layer=layer
        )

    def delete_layer_by_name(self, name):
        self.__structure.delete_layer_by_name(name)

    def set_source_layer_by_name(self, name):
        self.__structure.set_source_layer_by_name(name)

    def set_probe_layer_by_name(self, name):
        self.__target_layer_index = self.__structure.get_layer_index_by_name(name)

    def delete_layer_by_index(self, index):
        self.__structure.delete_layer_by_index(name)

    def set_source_layer_by_index(self, index):
        self.__structure.set_source_layer_by_index(index)

    def set_probe_layer_by_index(self, index):
        self.__target_layer_index = index

    def set_probe_layer_z_coordinate(self, target_z):
        self.__target_z = target_z

    def set_options(self, options):
        self.__options = options

    def output_sys_info(self):
        '''
        Print structure details, simulation options and settings
        '''
        pass

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

    def set_num_of_G(self, nG):
        self.__structure.init_structure(nG=nG, truncation=self.__options.truncation)

    def get_num_of_G(self):
        return self.__structure.get_nG()

    def get_G(self):
        return self.__structure.get_G()

    def get_phi(self):
        if self.__phi is None:
            raise Exception("Phi is None, Please integrate firstly!")
        return self.__phi

    def get_optical_parameters(self, omege_index, position):
        '''
        Get the optical parameters retrievaled by RCWA
        '''
        MICRON = 1e6 # micro meter unit
        position = MICRON*position # convert the SI unit (m) into the micro meter (um)
        if omege_index < 0:
            raise Exception("Omega index is smaller than zero!")
        if omege_index >= len(self.__omega_list) is None:
            raise Exception("Omega index of of range!")

        # initiate the RCWA for the given omega
        rcwa = calculation.RCWA(self.__structure, verbose=self.verbose)

        layer_index = 0
        offset = 0
        layers = self.__structure.get_layers()
        thickness_list = self.__structure.get_thickness_list()
        n_layer = len(layers)
        for i in range(n_layer):
            if (position[2] > offset) and (position[2] <= offset + thickness_list[i]):
                layer_index = i
                break
            offset += thickness_list[i]

        # only one layer exist
        if (layer_index == 0) and (position[2] > offset):
            layer_index = num_layer - 1

        Gx, Gy = self.__structure.get_G()
        nG = self.__structure.get_nG()
        Gx_r, Gx_l = np.meshgrid(Gx, Gx)
        Gy_r, Gy_l = np.meshgrid(Gy, Gy)

        Gx_mat = Gx_l - Gx_r
        Gy_mat = Gy_l - Gy_r
        r1, r2, r3, r4 = 0, nG, nG, 2*nG
        pos = (nG - 1)/2
        if (self.__options.truncation == "circular") and (self.__structure.get_dimension() == "two"):
            pos = 0
        
        x, y, z = position
        phase = np.exp(-1j*(Gx_mat[pos]*x + Gy_mat[0]*y))

        # TODO re-construct the real-space distribution of the optical parameters based on results from the Fourier transformations

    def get_phi_at_kx_ky(self, omega_index, kx, ky, target_layer_index, target_z):
        if (omega_index >= len(self.__omega_list)) or (omega_index < 0):
            raise Exception("Omega_index: Out of range! [{0} - {1}]".format(0, len(self.__omega_list)-1))

        rcwa = calculation.RCWA(self.__structure, verbose=self.verbose)
        poyntingflux = rcwa.get_poynting_flux(
            kx=kx, ky=ky, 
            omega_index=omega_index, 
            target_layer_index=target_layer_index, 
            target_z=target_z, 
            polarization=self.__options.polarization
        )
        c0 = 299792458.0 # the light speed in the vacuum
        omega = self.__omega_list[omega_index]
        phi = omega/c_0/(2*np.pi)*poyntingflux
        return phi

    def set_kx_integral(self, kpoints):

    def set_kx_integral_sym(self, points):

    def set_ky_integral(self, kpoints):

    def set_ky_integral_sym(self, points, end=0):

    def integrate_ky_ky(self, n_core=1):


        if n_core > 1:
            from multiprocessing import Pool
            from multiprocessing.dummy import Pool as ThreadPool


class SimulationPlanar(Simulation):
    pass

class SimulationGrating(Simulation):
    pass

class SimulationPattern(Simulation):
    pass