import numpy as np 

def mesh_grid(vL, vR):
    qL, qR = np.meshgrid(vL, vR)
    return qL, qR

def get_S_matrices(start_layer, N, num_of_layer, M_matrices, F_matrices, S_matrices, direction):
    r1 =0
    r2 = 2*N - 1
    r3 = 3*N
    r4 = 4*N - 1

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

def get_grand_imaginary_matrices(grand_imaginary_matrices, im_eps_xx, im_eps_xy, im_eps_yx, im_eps_yy, im_eps_zz, num_of_layer, N):

def get_E_matrices(E_matrices, eps_xx, eps_xy, eps_yx, eps_yy, num_of_layer, N):

def poynting_flux(omega, thickness_list, kx, ky, E_matrices, grand_imaginary_matrices, eps_zz_inv, Gx_mat, Gy_mat, source_list, target_layer, N, polar, z):
