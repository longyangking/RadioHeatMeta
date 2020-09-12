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
import radioheatmeta.material as Material

class FileLoader:
    def __init__(self, filename, verbose=False):
        self.filename = filename
        self.__omega_list = None
        self.__epsilon_list = None
        self.__mu_list = None
        self.__num_of_omega = 0
        self.__pre_set = False
        self.verbose = verbose

    def load_epsilon(self, filename):
        # Check type
        epsilon_type = 0
        with open(filename, 'r') as f:
            line = f.readline()
            epsilon_type = self.__check_type(line, filename)

        # Read data
        with open(filename, 'r') as f:
            line = f.readline()
            while line:
                line_str = line.strip().split()
                line_data = [float(x) for x in line_str]
                omega = line_data[0]
                self.__omega_list.append(omega)

                # Load epsilon data based on its type
                if epsilon_type == "scalar":
                    eps_re, eps_im = line_data[1:]
                    epsilon_vals = eps_re - 1j*eps_im
                    
                elif epsilon_type == "diagonal":
                    eps_re_xx, eps_im_xx, eps_re_yy, eps_im_yy, eps_re_zz, eps_im_zz = line_data[1:]
                    epsilon_vals = np.array([eps_re_xx - 1j*eps_im_xx, eps_re_yy - 1j*eps_im_yy, eps_re_zz - 1j*eps_im_zz], dtype=complex)

                else: # Tensor type
                    eps_re_xx, eps_im_xx, eps_re_yy, eps_im_yy, eps_re_zz, eps_im_zz, eps_re_xy, eps_im_xy, eps_re_xz, eps_im_xz, eps_re_yz, eps_im_yz = line_data[1:]
                    eps_xx = eps_re_xx - 1j*eps_im_xx
                    eps_yy = eps_re_yy - 1j*eps_im_yy
                    eps_zz = eps_re_zz - 1j*eps_im_zz
                    eps_xy = eps_re_xy - 1j*eps_im_xy
                    eps_xz = eps_re_xz - 1j*eps_im_xz
                    eps_yz = eps_re_yz - 1j*eps_im_yz
                    epsilon_vals = np.array([[eps_xx, eps_xy, eps_xz], [np.conj(eps_xy), eps_yy, eps_yz], [np.conj(eps_xz), np.conj(eps_yz), eps_zz]], dtype=complex)

                epsilon = Material.Epsilon(epsilon_vals, epsilon_type)
                self.__epsilon_list.append(epsilon)  

                line = f.readline()
    
    def load_mu(self, filename):
        # Check type
        mu_type = 0
        with open(filename, 'r') as f:
            line = f.readline()
            mu_type = self.__check_type(line, filename)

        # Read data
        with open(filename, 'r') as f:
            line = f.readline()
            while line:
                line_str = line.strip().split()
                line_data = [float(x) for x in line_str]
                omega = line_data[0]
                self.__omega_list.append(omega)

                # Load epsilon data based on its type
                if mu_type == "scalar":
                    mu_re, mu_im = line_data[1:]
                    mu_vals = mu_re - 1j*mu_im
                    
                elif epsilon_type == "diagonal":
                    mu_re_xx, mu_im_xx, mu_re_yy, mu_im_yy, mu_re_zz, mu_im_zz = line_data[1:]
                    mu_vals = np.array([mu_re_xx - 1j*mu_im_xx, mu_re_yy - 1j*mu_im_yy, mu_re_zz - 1j*mu_im_zz], dtype=complex)

                else: # Tensor type
                    mu_re_xx, mu_im_xx, mu_re_yy, mu_im_yy, mu_re_zz, mu_im_zz, mu_re_xy, mu_im_xy, mu_re_xz, mu_im_xz, mu_re_yz, mu_im_yz = line_data[1:]
                    mu_xx = mu_re_xx - 1j*mu_im_xx
                    mu_yy = mu_re_yy - 1j*mu_im_yy
                    mu_zz = mu_re_zz - 1j*mu_im_zz
                    mu_xy = mu_re_xy - 1j*mu_im_xy
                    mu_xz = mu_re_xz - 1j*mu_im_xz
                    mu_yz = mu_re_yz - 1j*mu_im_yz
                    mu_vals = np.array([[mu_xx, mu_xy, mu_xz], [np.conj(mu_xy), mu_yy, mu_yz], [np.conj(mu_xz), np.conj(mu_yz), mu_zz]], dtype=complex)

                mu = Material.Mu(mu_vals, mu_type)
                self.__mu_list.append(mu)  

                line = f.readline()

    def __check_type(self, line, filename):
        line_str = line.strip().split()
        #line_data = [float(x) for x in line_str]
        num_of_space = len(line_str) - 1 # First string is the frequency
        data_type = None
        if num_of_space == 2:
            data_type = "scalar"
        elif num_of_space == 6:
            data_type = "diagonal"
        elif num_of_space == 10:
            data_type = "tensor"
        else:
            raise Exception("Wrong input type in {filename}: should be of 2, 6, 10 numbers for each frequency!".format(filename=filename))

        return data_type

    def get_omega_list(self):
        return self.__omega_list

    def get_epsilon_list(self):
        return self.__epsilon_list

    def get_mu_list(self):
        return self.__mu_list

    def get_num_of_omega(self):
        return len(self.__omega_list)
