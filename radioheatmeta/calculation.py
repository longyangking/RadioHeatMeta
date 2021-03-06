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
import scipy as sp 

class RCWA:
    def __init__(self, structure, verbose=False):
        self.__structure = structure
        self.verbose = verbose


        self.__S_matrices = None

    def init_RCWA(self, omega_index):


    def get_grand_imaginary_matrices(self, omega_index):
        '''
        Get the imaginary part of optical parameters
        '''
        N = self.get_nG()
        Gx, Gy = self.__G
        num_of_layer = len(self.__layer_list)

        grand_imaginary_matrices = list()
        for layer in self.__layer_list:
            [epsxx, epsxy, epsxz, \
            epsyx, epsyy, epsyz, \
            epszx, epszy,  epszz, \
            muxx, muxy, muxz, \
            muyx, muyy, muyz, \
            muzx, muzy, muzz] = layer.get_optical_matrices(Gx, Gy, omega_index):

            grand_imaginary_epsilon_matrix = np.zeros((3*N, 3*N), dtype=complex)
            grand_imaginary_epsilon_matrix[0:N, 0:N] = np.imag(epsxx)
            grand_imaginary_epsilon_matrix[0:N, N:2*N] = np.imag(epsxy)
            grand_imaginary_epsilon_matrix[0:N, 2*N:3*N] = np.imag(epsxz)
            grand_imaginary_epsilon_matrix[N:2*N, 0:N] = im_eps_yx[i]
            grand_imaginary_epsilon_matrix[N:2*N, N:2*N] = im_eps_yy[i]
            grand_imaginary_epsilon_matrix[2*N:3*N, 2*N:3*N] = im_eps_zz[i]

            grand_imaginary_mu_matrix = np.zeros((3*N, 3*N), dtype=complex)
            grand_imaginary_mu_matrix[0:N, 0:N] = np.imag(muxx)
            grand_imaginary_mu_matrix[0:N, N:2*N] = np.imag(epsxy)
            grand_imaginary_mu_matrix[0:N, 2*N:3*N] = np.imag(epsxz)
            grand_imaginary_mu_matrix[N:2*N, 0:N] = im_eps_yx[i]
            grand_imaginary_mu_matrix[N:2*N, N:2*N] = im_eps_yy[i]
            grand_imaginary_mu_matrix[2*N:3*N, 2*N:3*N] = im_eps_zz[i]

            grand_imaginary_matrices.append([grand_imaginary_epsilon_matrix, grand_imaginary_mu_matrix])
            
        return grand_imaginary_matrices

    def get_E_matrices(self, omega_index):
        '''
        Get E Matrices based on the reciprocal lattice matrices
        '''
        nG = self.get_nG()

        # E_matrices = np.zeros((num_of_layer, 2*N, 2*N))
        # for i in range(num_of_layer):
        #     E_matrix = np.concatenate((
        #         np.concatenate((eps_yy[i], -eps_yx[i]), axis=1),
        #         np.concatenate((-eps_xy[i], eps_xx[i]), axis=1)
        #     ), axis=0)
        #     E_matrices[i] = E_matrix
        # return E_matrices

    def get_M_matrices(self, omega_index):
        '''
        Get M Matrices based on the reciprocal lattice matrices
        '''
        nG = self.get_nG()


    def get_S_matrices(self):
        '''
        Get the scattering matrices S
        '''
        r1 =0
        r2 = 2*N
        r3 = 3*N
        r4 = 4*N

        if (direction == "all") or (direction == "down"):
            for i in range(start_layer,0,-1):
                I = np.linalg.solve(M_matrices[i], M_matrices[i-1])
                
                left_top = I[r3:r4, r3:r4]
                right_top = I[r3:r4, r1:r2]
                left_bottom = I[r1:r2, r3:r4]
                right_bottom = I[r1:r2, r1:r2]

                S_matrices[i-1][r1:r2, r1:r2] = np.linalg.solve(
                    left_top - F_matrices[i].dot(S_matrices[i][r1:r2, r3:r4].dot(left_bottom)),
                    F_matrices[i]
                ).dot(S_matrices[i][r1:r2, r1:r2])

                test = left_top - F_matrices[i].dot(S_matrices[i][r1:r2, r3:r4].dot(left_bottom))

                S_matrices[i-1]

    def get_poynting_flux(self, kx, ky, omega_index, target_layer_index, target_z, polarization="both"):
        kx = kx*omega
        ky = ky*omega
        r1 = 0
        r2 = 2*N-1
        r3 = 2*N
        r4 = 4*N-1
        one_padding_4N = np.diag(np.ones(4*N), dtype=complex)
        one_padding_2N = np.diag(np.ones(2*N), dtype=complex)
        one_padding_1N = np.diag(np.ones(N), dtype=complex)
        zero_padding_2N = np.zeros((2*N, 2*N), dtype=complex)
        zero_padding_4N = np.zeros((4*N, 4*N), dtype=complex)
        num_of_layer = len(thickness_list)

        kx_mat = np.diag(kx + Gx_mat)
        ky_mat = np.diag(ky + Gy_mat)

        T_matrices = list()
        M_matrices = list()
        Eigen_val_matrices = list()
        Eigen_vec_matrices = list()
        F_matrices = list()
        coeff_of_A = one_padding_2N
        coeff_of_B = one_padding_2N

        k_matrix = np.concatenate(
            (
            np.concatenate((kx_mat.dot(kx_mat), kx_mat.dot(ky_mat)), axis=1),
            np.concatenate((ky_mat.dot(kx_mat), ky_mat.dot(ky_mat)), axis=1)   
            ),
            axis=0
        )

        for i in range(num_of_layer):
            vertical_align = np.concatenate((ky_mat, -kx_mat), axis=0)
            horizontal_align = np.concatenate((ky_mat, kx_mat), axis=1)


        return flux