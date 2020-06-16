import numpy as np 

def get_G_matrices(n_G, reciprocal_lattice, dimension, truncation):
    if (n_G <=0) and (dimension != "no"):
        raise Exception("Require number of G more than 1!")
    if dimension not in ["no", "one", "two"]:
        raise Exception("Unknown dimension type: {dimension}. Only support: no, one or two.".format(dimension=dimension))
    if truncation not in ["circular", "parallelogramic"]:
        raise Exception("Unknown truncation type: {truncation}. Only support: circular, parallelogramic".format(truncation=truncation))

    if (dimension == "no"):
        Gx_mat, Gy_mat, n_G = __get_G_parallelogramic(n_G=1, reciprocal_lattice, dimension)
    elif (dimension == "one"):

    else:
        if (truncation == "circular"):

        else:

    return Gx_mat, Gy_mat, n_G

def __get_G_parallelogramic(n_G, reciprocal_lattice, dimension):
    if (dimension == "one"): # 1D case
        m = int(n_G/2)
        n_G = 2*m + 1
        Gx_list, Gy_list = np.zeros(n_G), np.zeros(1)
        for i in range(-M, M+1):
            Gx_list(i+M) = i
        Gx_mat, Gy_mat = np.meshgrid(Gx_list, Gy_list)
    else: # 2D case
        n_root = int(np.sqrt(n_G))
        if (n_root % 2 ==0) and (n_root > 0):
            n_root -= 1
        m = int(n_root / 2)

        G_list = np.zeros(n_root)
        for i in range(-M, M+1):
            G_list(i+M) = i
        Gx_mat, Gy_mat = np.meshgrid(G_list, G_list)
        n_G = np.power(n_root, 2)

    Gx_mat_temp = reciprocal_lattice.bx[0] * Gx_mat + reciprocal_lattice.by[0] * Gy_mat
    Gy_mat = reciprocal_lattice.bx[1] * Gx_mat + reciprocal_lattice.by[1] * Gy_mat
    Gx_mat = Gx_mat_temp
    return Gx_mat, Gy_mat, n_G
