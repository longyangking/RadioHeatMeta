import numpy as np 
from enum import Enum
import geometry as geom

class Dimension_type(Enum):
    No = 0
    One = 1
    Two = 2

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

        if epsilon_type not in ["scalar", "diagonal", "tensor"]:
            raise Exception("Wrong epsilon type: [{epsilon_type}]".format(epsilon_type=epsilon_type)) 
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

        if mur_type not in ["scalar", "diagonal", "tensor"]:
            raise Exception("Wrong mur type: [{mur_type}]".format(mur_type=mur_type)) 
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

class Pattern:
    def __init__(self, args1=[0,0], args2=[0,0], pattern_type="rectangle", area=0, edge_list=None, parent=-1, angle=0):
        self.args1 = args1
        self.args2 = args2
        self.edge_list = edge_list

        if pattern_type not in ["grating", "rectangle", "circle", "ellipse", "polygon"]:
            raise Exception("Wrong pattern type: [{patter_type}]".format(pattern_type=pattern_type))
        self.pattern_type = pattern_type

        self.area = area
        self.parent = parent
        self.angle = angle

class Layer:
    def __init__(self, name, background_material, thickness):
        self.__name = name 
        self.__thickness = thickness
        self.__has_tensor = False

        self.__is_source = False
        self.__background_material = background_material

        self.__material_list = list()
        self.__pattern_list = list()
    
    @staticmethod
    def copy(name):
        # TODO copy material
        pass

    def set_background(self, material):
        self.__background_material = material
        epsilon_type, mur_type = material.material_type
        if (epsilon_type == "tensor") or (mur_type == "tensor"):
            self.__has_tensor = True

    def get_background(self):
        return self.__background_material

    def set_thickness(self, thickness):
        self.__thickness = thickness

    def get_thickness(self):
        return self.__thickness

    def set_is_source(self, status=True):
        self.__is_source = status

    def check_is_source(self):
        return self.__is_source

    def contain_tensor(self, status=False):
        self.__has_tensor = status

    def has_tensor(self):
        return self.__has_tensor

    def has_material(self, material):
        if (self.__background_material == material):
            return True
        for material in self.__material

    def get_material_by_name(self, name):
        if len(self.__material_list) == 0:
            raise Exception("Zero length for the material list for the layer \"{name}\"".format(name=self.__name))
        for i in range(len(self.__material_list)):
            material = self.__material_list[i]
            if material.name == name:
                return i
        return None

    def get_num_of_material(self):
        return len(self.__material_list)

    def get_name(self):
        return self.__name

    def get_materials_begin(self):
        pass

    def get_materials_end(self):
        pass

    def get_material_list(self):
        return self.__material_list

    def get_patterns_begin(self):
        pass

    def get_patterns_end(self):
        pass

    def get_pattern_list(self):
        return self.__pattern_list

    def add_rectangle_pattern(self, material, args1, args2, angle):
        '''
            args1: the position of centers (x,y)
            angle: the rotated angle with respect to x axis
            args2: the widths in x and y directions
        '''
        self.__material_list.append(material)
        epsilon_type, mur_type = material.material_type
        if (epsilon_type == "tensor") or (mur_type == "tensor"):
            self.__has_tensor = True
        pattern = Pattern(
            args1=args1, 
            args2=args2,
            angle=angle, 
            type="rectangle",
            area= geom.get_rectangle_area(args2)
        )
        self.__pattern_list.append(pattern)

    def add_circle_pattern(self, material, args, radius):
        '''
        args1: the position of centers (x,y)
        radius: the radius of the circle
        '''
        self.__material_list.append(material)
        epsilon_type, mur_type = material.material_type
        if (epsilon_type == "tensor") or (mur_type == "tensor"):
            self.__has_tensor = True
        pattern = Pattern(
            args1 = [args[0], radius],
            args2 = [args[1], radius],
            pattern_type="circle",
            area = geom.get_circle_area(radius)
        )
        self.__pattern_list.append(pattern)

    def add_ellipse_pattern(self, material, args1, angle, args2):
        '''
        args1: the position of centers (x,y)
        angle: the rotated angle with respect to x axis
        args2: the halfwidths in x and y directions
        '''
        self.__material_list.append(Material)
        epsilon_type, mur_type = material.material_type
        if (epsilon_type == "tensor") or (mur_type == "tensor"):
            self.__has_tensor = True
        pattern = Pattern(
            args1=args1,
            args2=args2,
            angle=angle,
            pattern_type="ellipse",
            area=geom.get_ellipse_area(args2)
        )
        self.__pattern_list.append(pattern)

    def add_polygon_pattern(self, material, args1, angle, edge_points):
        '''
        args1: the position of centers (x,y)
        angle: the rotated angle with respect to x axis
        edgePoints: the points of the vertices in counter clockwise order
        '''
        self.__material_list.append(material)
        epsilon_type, mur_type = material.material_type
        if (epsilon_type == "tensor") or (mur_type == "tensor"):
            self.__has_tensor = True
        pattern = Pattern(
            args1=args1,
            angle=angle,
            pattern_type="polygon",
            edge_points=edge_points,
            area=geom.get_polygon_area(edge_list)
        )
        self.__pattern_list.append(pattern)

    def add_grating_pattern(self, material, center, width):
        
        self.__material_list.append(material)
        epsilon_type, mur_type = material.material_type
        if (epsilon_type == "tensor") or (mur_type == "tensor"):
            self.__has_tensor = True
        pattern = Pattern(
            args1=[center, width],
            args2=[0,0],
            pattern_type="grating",
            area=geom.get_grating_area(width)
        )
        self.__pattern_list.append(pattern)

    def __is_contain_in_geometry(self, pattern1, pattern2):
        center2 = np.zeros(2)
        patten_type = pattern1.pattern_type
        


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