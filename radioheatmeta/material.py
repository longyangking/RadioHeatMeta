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

#  Define several operations for type conversions or Fourier transformations of material data 

import numpy as np 
from radioheatmeta.fileloader import FileLoader

class Epsilon:
    def __init__(self, epsilon_val, epsilon_type="scalar"):
        self.__epsilon_val = np.array(epsilon_val, dtype=complex)

        if epsilon_type not in ["scalar", "diagonal", "tensor"]:
            raise Exception("Wrong epsilon type: [{epsilon_type}]".format(epsilon_type=epsilon_type)) 
        self.__epsilon_type = epsilon_type

        self.__init_vals()

    def __init_vals(self):
        '''
        Transform the all data type into the tensor
        '''
        if self.__epsilon_type == "scalar":
            epsilon_val = np.diag(self.__epsilon_val * np.ones(3,dtype=complex), dtype=complex)
        elif self.__epsilon_type == "diagonal":
            epsilon_val = np.diag(self.__epsilon_val, dtype=complex)
        
        self.__epsilon_val = epsilon_val

    @property
    def epsilon_val(self):
        return self.__epsilon_val

    @property
    def epsilon_type(self):
        return self.__epsilon_type

class Mu:
    def __init__(self, mu_val=1, mu_type="scalar"):
        self.__mu_val = np.array(mu_val, dtype=complex)

        if mu_type not in ["scalar", "diagonal", "tensor"]:
            raise Exception("Wrong mu type: [{mu_type}]".format(mu_type=mu_type)) 
        self.__mu_type = mu_type

        self.__init_vals()

    def __init_vals(self):
        '''
        Transform the all data type into the tensor
        '''
        if self.__mu_type == "scalar":
            mu_val = np.diag(self.__mu_val * np.ones(3,dtype=complex), dtype=complex)
        elif self.__mu_type == "diagonal":
            mu_val = np.diag(self.__mu_val, dtype=complex)
        
        self.__mu_val = mu_val

    @property
    def mu_val(self):
        return self.__mu_val

    @property
    def mu_type(self):
        return self.__mu_type

class Material:
    '''
    Material class for multilayer systems
    '''
    def __init__(self, name, omega_list, epsilon_list, mu_list):
        self.__name = name
        self.__omega_list = omega_list
        self.__epsilon_list = epsilon_list
        self.__mu_list = mu_list

        if len(self.__epsilon_list) == 0:
            raise Exception("Zero vals in the epsilon list for the material \"{name}\"".format(name=self.__name))
        if len(self.__mu_list) == 0:
            raise Exception("Zero vals in the mu list for the material \"{name}\"".format(name=self.__name))

    def __init__(self, name, omega_list, epsilon_vals, mu_vals, epsilon_type="scalar", mu_type="scalar"):
        self.__name = name
        self.__omega_list = omega_list
        self.__epsilon_list = list()
        self.__mu_list = list()

        for i in range(len(self.__omega_list)):
            epsilon = Epsilon(epsilon_val=epsilon_vals[i], epsilon_type=epsilon_type)
            mu = Mu(mu_val=mu_vals[i], mu_type=mu_type)
            self.__epsilon_list.append(epsilon)
            self.__mu_list.append(mu)
    
    @property
    def name(self):
        '''
        Get the name of the material
        '''
        return self.__name

    def get_omega_list(self):
        return self.__omega_list

    def get_material_type(self):
        if len(self.__epsilon_list) == 0:
            raise Exception("Zero vals in the epsilon list for the material \"{name}\"".format(name=self.__name))
        if len(self.__mu_list) == 0:
            raise Exception("Zero vals in the mu list for the material \"{name}\"".format(name=self.__name))

        epsilon = self.__epsilon_list[0]
        mu = self.__mu_list[0]
        return epsilon.epsilon_type, mu.mu_type

    def get_epsilon_at_index(self, index):
        if index < 0:
            raise Exception("Index: smaller than zero for the material \"{name}\"".format(name=self.__name))
        if len(self.__epsilon_list) <= index:
            raise Exception("Index: larger than the length of the epsilon list for the material \"{name}\"".format(name=self.__name))

        return self.__epsilon_list[index].epsilon_val

    def get_mu_at_index(self, index):
        if index < 0:
            raise Exception("Index: smaller than zero for the material \"{name}\"".format(name=self.__name))
        if len(self.__mu_list) <= index:
            raise Exception("Index: larger than the length of the epsilon list for the material \"{name}\"".format(name=self.__name))

        return self.__mu_list[index].mu_val

    def get_epsilon_list(self):
        return self.__epsilon_list

    def get_mu_list(self):
        return self.__mu_list

## TODO
def load_material(name, filename, verbose=False):
    fileloader = FileLoader(filename, verbose=verbose)
    material = Material(
        name=name,
        omega_list=fileloader.get_omega_list(),
        epsilon_list=fileloader.get_epsilon_list(),
        mur_list=fileloader.get_mur_list()
    )
    if self.verbose:
        print("import material:[{name}] into the simulation".format(name=name))

    self.__material_map.append([name, material])
    self.__structure.add_material(material)