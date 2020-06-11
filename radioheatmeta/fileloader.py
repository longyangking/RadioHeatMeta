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
from system import Epsilon, Mur

class FileLoader:
    def __init__(self, verbose=False):
        self.__omega_list = None
        self.__epsilon_list = None
        self.__mur_list = None
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

                epsilon = Epsilon(epsilon_vals, epsilon_type)
                self.__epsilon_list.append(epsilon)  

                line = f.readline()
    
    def load_mur(self, filename):
        # Check type
        mur_type = 0
        with open(filename, 'r') as f:
            line = f.readline()
            mur_type = self.__check_type(line, filename)

        # Read data
        with open(filename, 'r') as f:
            line = f.readline()
            while line:
                line_str = line.strip().split()
                line_data = [float(x) for x in line_str]
                omega = line_data[0]
                self.__omega_list.append(omega)

                # Load epsilon data based on its type
                if mur_type == "scalar":
                    mur_re, mur_im = line_data[1:]
                    mur_vals = mur_re - 1j*mur_im
                    
                elif epsilon_type == "diagonal":
                    mur_re_xx, mur_im_xx, mur_re_yy, mur_im_yy, mur_re_zz, mur_im_zz = line_data[1:]
                    mur_vals = np.array([mur_re_xx - 1j*mur_im_xx, mur_re_yy - 1j*mur_im_yy, mur_re_zz - 1j*mur_im_zz], dtype=complex)

                else: # Tensor type
                    mur_re_xx, mur_im_xx, mur_re_yy, mur_im_yy, mur_re_zz, mur_im_zz, mur_re_xy, mur_im_xy, mur_re_xz, mur_im_xz, mur_re_yz, mur_im_yz = line_data[1:]
                    mur_xx = mur_re_xx - 1j*mur_im_xx
                    mur_yy = mur_re_yy - 1j*mur_im_yy
                    mur_zz = mur_re_zz - 1j*mur_im_zz
                    mur_xy = mur_re_xy - 1j*mur_im_xy
                    mur_xz = mur_re_xz - 1j*mur_im_xz
                    mur_yz = mur_re_yz - 1j*mur_im_yz
                    mur_vals = np.array([[mur_xx, mur_xy, mur_xz], [np.conj(mur_xy), mur_yy, mur_yz], [np.conj(mur_xz), np.conj(mur_yz), mur_zz]], dtype=complex)

                mur = Mur(mur_vals, mur_type)
                self.__mur_list.append(mur)  

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

    def get_mur_list(self):
        return self.__mur_list

    def get_num_of_omega(self):
        return len(self.__omega_list)
