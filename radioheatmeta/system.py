import numpy as np 
from enum import Enum

class Dimension_type(Enum):
    No = 0
    One = 1
    Two = 2

class Pattern_type(Enum):
    Grating = 0
    Rectangle = 1
    Circle = 2
    Ellipse = 3
    Polygon = 4

class Material_type(Enum):
    Scalar = 0
    Diagonal = 1
    Tensor = 2

class Polarization_type(Enum):
    TE = 0
    TM = 1
    Both = 2

class Truncation_type(Enum):
    Circular = 0
    Parallelogramic = 1

class Lattice:
    def __init__(self, bx=[0,0], by=[0,0], angle=90, area=0):
        self.__bx = np.array(bx)
        self.__by = np.array(by)
        self.__angle = angle
        self.__area = area

    @property
    def bx(self):
        return self.__bx
    
    @bx.setter
    def bx(self, bx):
        self.__bx = bx

    @property
    def by(self):
        return self.__by

    @by.setter
    def by(self):
        self.__by = by

    @property
    def angle(self):
        return self.__angle

    @angle.setter
    def angle(self, angle):
        self.__angle = angle

    @property
    def area(self):
        return self.__area

    @area.setter
    def area(self, area):
        self.__area = area

class Epsilon:
    def __init__(self, epsilon_vals, epsilon_type="scalar"):
        self.__epsilon_vals = np.array(epsilon_vals)
        self.__epsilon_type = epsilon_type

    @property
    def epsilon_vals(self):
        return self.__epsilon_vals

    @epsilonvals.setter
    def epsilon_vals(self, epsilon_vals):
        self.__epsilon_vals = epsilon_vals

    @property
    def epsilon_type(self):
        return self.__epsilon_type
    
    @epsilon_type.setter
    def epsilon_type(self, epsilon_type):
        self.__epsilon_type = epsilon_type

class Mur:
    def __init__(self, mur_vals, mur_type="scalar"):
        self.__mur_vals = mur_vals 
        self.__mur_type = mur_type

    @property
    def mur_vals(self):
        return self.__murvals

    @mur_vals.setter
    def mur_vals(self, epsilonvals):
        self.__murvals = mur_vals

    @property
    def mur_type(self):
        return self.__mur_type

    @mur_type.setter
    def mur_type(self, mur_type):
        self.__mur_type = mur_type

class Material:
    def __init__(self, name, omega_list, epsilon_list, mur_list):
        self.__name = name
        self.__omega_list = omega_list
        self.__epsilon_list = epsilon_list
        self.__mur_list = mur_list

        if len(self.__epsilon_list) == 0:
            raise Exception("Zero vals in the epsilon list for the material \"{name}\"".format(name=self.__name))
        if len(self.__mur_list) == 0:
            raise Exception("Zero vals in the mur list for the material \"{name}\"".format(name=self.__name))
    
    @property
    def name(self):
        '''Get the name of the material

        Args:
            No

        Returns:
            Name (string)
        '''
        return self.__name

    @property
    def material_type(self):
        if len(self.__epsilon_list) == 0:
            raise Exception("Zero vals in the epsilon list for the material \"{name}\"".format(name=self.__name))
        if len(self.__mur_list) == 0:
            raise Exception("Zero vals in the mur list for the material \"{name}\"".format(name=self.__name))

        epsilon = self.__epsilon_list[0]
        mur = self.__mur_list[0]
        return epsilon.type, mur.type

    def get_epsilon_at_index(self, index):
        if index < 0:
            raise Exception("Index: smaller than zero for the material \"{name}\"".format(name=self.__name))
        if len(self.__epsilon_list) <= index:
            raise Exception("Index: larger than the length of the epsilon list for the material \"{name}\"".format(name=self.__name))

        return self.__epsilon_list[index].epsilon_vals

    @property
    def omega_list(self):
        return self.__omega_list

    @omegas.setter
    def omega_list(self, omega_list):
        self.__omega_list = omega_list

    @property
    def epsilon_list(self):
        return self.__epsilon_list

    @epsilons.setter
    def epsilon_list(self, epsilon_list):
        # TODO Check before setting
        self.__epsilon_list = epsilon_list

class Layer:
    def __init__(self, name, material, thickness):
        self.__name = name 
        self.__material = material
        self.__thickness = thickness
        self.__has_tensor = False
    
    @staticmethod
    def copy(name):
        pass

    def set_background(self, material):
        pass
        
    def get_background(self):
        pass

    def set_thickness(self, thickness):
        pass

    def get_thickness(self):
        pass

    def set_is_source(self):
        pass

    def check_is_source(self):
        pass

    def contain_tensor(self, status):
        pass

    def has_tensor(self):
        pass

    def has_material(self, material):
        pass

    def get_material_by_name(self, name):
        pass

    def get_num_of_material(self):
        pass

    def get_name(self):
        pass

    def get_materials_begin(self):

    def get_materials_end(self):

    def get_patterns_begin(self):

    def get_patterns_end(self):

    def add_rectangle_pattern(self, material, args1, args2, angle):

    def add_circle_pattern(self, material, args, radius):

    def add_ellipse_pattern(self, material, args1, angle, args2):

    def add_polygon_pattern(self, material, args1, angle, edge_points):

    def add_grating_pattern(self, material, center, width):

    def get_geometry_containment_relation(self):


class Structure:
    def __init__(self):

    def __init__(self, structure):

    def add_material(self, material):

    def add_layer(self, layer):

    def set_lattice(self, lattice):

    def delete_layer_by_name(self, name):

    def delete_layer_by_layer(self, layer):

    def get_layer_by_index(self, index):

    def get_layer_by_name(self, name):

    def get_num_of_layer(self):

    def get_thickness_list(self):

    def get_layers_begin(self):

    def get_layers_end(self):

    def delete_layer(self, it):

    def reorganize_layers(self):