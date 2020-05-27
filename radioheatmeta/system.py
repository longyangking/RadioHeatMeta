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
    def __init__(self, epsilonvals):
        self.__epsilonvals = epsilonvals

    @property
    def epsilonvals(self):
        return self.__epsilonvals

    @epsilonvals.setter
    def epsilonvals(self, epsilonvals):
        self.__epsilonvals = epsilonvals

class Mur:
    def __init__(self, murvals):
        self.__murvals = murvals 

    @property
    def murvals(self):
        return self.__murvals

    @murvals.setter
    def murvals(self, epsilonvals):
        self.__murvals = murvals

class Material:
    def __init__(self, name, omegas, epsilons):
        self.name = name
    
    @property
    def name(self):
        return self.name

    @property
    def material_type(self):
        return self.material_type

    def get_epsilon_at_index(self, index):
        pass

    @property
    def omegas(self):
        return self.omegas

    @omegas.setter
    def omegas(self, omegas):
        self.omegas = omegas

    @property
    def epsilons(self):
        return self.epsilons

    @epsilons.setter
    def epsilons(self, epsilons):
        # TODO Check before setting
        self.epsilons = epsilons
    
class Pattern:
    pass

class Layer:
    def __init__(self, name, material, thickness):
        self.name = name 
        self.material = material
        self.thickness = thickness
    
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

