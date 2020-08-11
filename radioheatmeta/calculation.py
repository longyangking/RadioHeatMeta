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

def sinc(x):
    if np.isscalar(x):
        if x == 0:
            return 1
        else:
            return np.sin(x)/x
    else:
        sin_x = np.sin(x)
        x = np.array(x)
        positions = np.where(x==0)
        x[positions] = 1
        result = sin_x / x
        result[positions] = 1
        return result

def jinc(x):
    if x == 0.0:
        return 0.5
    return sp.special.j1(x)/x

class RCWA:
    def __init__(self, structure):

    def __get_S_matrices(self):

    def __get_grand_imagnary_matrices(self):

    def __get_E_matrices(self):

    def get_poynting_flux(self, omega):

        return flux

def mesh_grid(vL, vR):
    qL, qR = np.meshgrid(vL, vR)
    return qL, qR

def get_S_matrices(start_layer, N, num_of_layer, M_matrices, F_matrices, S_matrices, direction):
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

def get_grand_imaginary_matrices(im_eps_xx, im_eps_xy, im_eps_yx, im_eps_yy, im_eps_zz, num_of_layer, N):
    '''
    Get the imaginary part of optical parameters 
    '''
    grand_imaginary_matrices = np.zeros((num_of_layer, 3*N, 3*N),dtype=complex)
    for i in range(num_of_layer):
        grand_imaginary_matrix = np.zeros((3*N, 3*N), dtype=complex)
        grand_imaginary_matrix[0:N, 0:N] = im_eps_xx[i]
        grand_imaginary_matrix[0:N, N:2*N] = im_eps_xy[i]
        grand_imaginary_matrix[N:2*N, 0:N] = im_eps_yx[i]
        grand_imaginary_matrix[N:2*N, N:2*N] = im_eps_yy[i]
        grand_imaginary_matrix[2*N:3*N, 2*N:3*N] = im_eps_zz[i]
        grand_imaginary_matrices[i] = grand_imaginary_matrix
    return grand_imaginary_matrices

def get_E_matrices(eps_xx, eps_xy, eps_yx, eps_yy, num_of_layer, N):
    '''
    Get Epsilon Matrices based on the reciprocal lattice matrices
    '''
    E_matrices = np.zeros((num_of_layer, 2*N, 2*N))
    for i in range(num_of_layer):
        E_matrix = np.concatenate((
            np.concatenate((eps_yy[i], -eps_yx[i]), axis=1),
            np.concatenate((-eps_xy[i], eps_xx[i]), axis=1)
        ), axis=0)
        E_matrices[i] = E_matrix
    return E_matrices

def poynting_flux(omega, thickness_list, kx, ky, E_matrices, grand_imaginary_matrices, eps_zz_inv, Gx_mat, Gy_mat, source_list, target_layer, N, polar, z):
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

        