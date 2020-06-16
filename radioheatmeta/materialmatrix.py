import numpy as np 
from system import Epsilon, Mur

def epsilon_scalar_2_diagonal(epsilon):
    if epsilon.epsilon_type != "scalar":
        raise Exception("Input epsilon is not scalar!")
    epsilon_vals = epsilon.epsilon_vals * np.ones(3,dtype=complex)
    return Epsilon(
        epsilon_vals=epsilon_vals,
        epsilon_type="diagonal"
    )

def mur_scalar_2_diagonal(mur):
    if mur.mur_type != "scalar":
        raise Exception("Input mur is not scalar!")
    mur_vals = mur.mur_vals * np.ones(3,dtype=complex)
    return Mur(
        mur_vals=mur_vals,
        mur_type="diagonal"
    )
    

def epsilon_diagonal_2_tensor(epsilon):
    if epsilon.epsilon_type != "diagonal":
        raise Exception("Input epsilon is not diagonal!")
    epsilon_vals = np.diag(epsilon.epsilon_vals)
    return Epsilon(
        epsilon_vals=epsilon_vals,
        epsilon_type="tensor"
    )

def mur_diagonal_2_tensor(mur):
    if mur.mur_type != "diagonal":
        raise Exception("Input mur is not diagonal!")
    mur_vals = np.diag(mur.mur_vals)
    return Mur(
        mur_vals=mur_vals,
        mur_type="tensor"
    )

def epsilon_scalar_2_tensor(epsilon):
    if epsilon.epsilon_type != "scalar":
        raise Exception("Input epsilon is not scalar!")
    epsilon_vals = np.diag(epsilon.epsilon_vals * np.ones(3,dtype=complex))
    return Epsilon(
        epsilon_vals=epsilon_vals,
        epsilon_type="tensor"
    )

def mur_scalar_2_tensor(mur):
    if mur.mur_type != "scalar":
        raise Exception("Input mur is not scalar!")
    mur_vals = np.diag(mur.mur_vals * np.ones(3,dtype=complex))
    return Mur(
        mur_vals=mur_vals,
        mur_type="tensor"
    )

def epsilon_to_tensor(epsilon):
    if epsilon.epsilon_type == "scalar":
        return epsilon_scalar_2_tensor(epsilon)
    elif epsilon.epsilon_type == "diagonal":
        return epsilon_diagonal_2_tensor(epsilon)

    return epsilon

def mur_to_tensor(mur):
    if mur.mur_type == "scalar":
        return mur_scalar_2_tensor(mur)
    elif mur.mur_type == "diagonal":
        return mur_diagonal_2_tensor(mur)

    return mur



