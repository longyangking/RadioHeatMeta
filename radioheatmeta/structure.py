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
from enum import Enum
import radioheatmeta.geometry as geometry
from radioheatmeta.lattice import Lattice
from radioheatmeta.material import Epsilon, Mu, Material

# class Polarization_type(Enum):
#     TE = 0
#     TM = 1
#     Both = 2

MICRON = 1e6

class Pattern:
    '''
    Define the geometrical pattern 
    '''
    def __init__(self, area, pattern_type="rectangle",  parent=-1):
        if pattern_type not in ["grating", "rectangle", "circle", "ellipse", "polygon"]:
            raise Exception("Wrong pattern type: [{patter_type}]".format(pattern_type=pattern_type))
        self.__pattern_type = pattern_type
        self.area = area
        self.__parent = parent

    def get_pattern_type(self):
        return self.__pattern_type

    def get_parent(self):
        return self.__parent

    def set_parent(self, parent):
        self.__parent = parent

class RectanglePatten(Pattern):
    def __init__(self, position, angle, widths, area):
        super().__init__(patter_type="rectangle", area=area)
        self.position = np.array(position)
        self.angle = angle
        self.widths = widths

class EllipsePattern(Pattern):
    def __init__(self, position, angle, halfwidths, area):
        super().__init__(patter_type="ellipse", area=area)
        self.position = np.array(position)
        self.angle = angle
        self.halfwidths = halfwidths

class PolygonPattern(Pattern):
    def __init__(self, position, angle, edge_list, area):
        super().__init__(pattern_type="polygon", area=area)
        self.position = np.array(position)
        self.angle = angle
        self.edge_list = edge_list

class CirclePattern(Pattern):
    def __init__(self, position, angle, area):
        super().__init__(pattern_type="circle", area=area)
        self.position = np.array(position)
        self.angle = angle

class GratingPattern(Pattern):
    '''
    Quasi 1D structure
    '''
    def __init__(self, center, width, area):
        super().__init__(pattern_type="grating", area=area)
        self.center = center
        self.width = width

class Layer:
    '''
    Define the information of each layer
    '''
    def __init__(self, 
        name, thickness, background_material,
        is_source=False, verbose=False
    ):
        self.__name = name 
        self.__thickness = thickness
        self.__has_tensor = False

        self.__is_source = is_source
        self.__background_material = background_material

        self.__material_list = list()
        self.__pattern_list = list()

        self.verbose = verbose
    
    @staticmethod
    def copy(name):
        # TODO copy material
        pass

    @property
    def name(self):
        return self.__name

    def set_background(self, material):
        self.__background_material = material
        epsilon_type, mu_type = material.get_material_type()
        if (epsilon_type == "tensor") or (mu_type == "tensor"):
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

    def contain_tensor(self, status=True):
        if not self.__has_tensor:
            self.__has_tensor = status

    def has_tensor(self):
        return self.__has_tensor

    def has_material(self, material):
        if (self.__background_material.name  == material.name):
            return True
        for imaterial in self.__material_list:
            if imaterial.name == material.name:
                return True
            
    def get_material_by_name(self, name):
        if len(self.__material_list) == 0:
            raise Exception("Zero length for the material list for the layer \"{name}\"".format(name=self.__name))
        for material in self.__material_list:
            if material.name == name:
                return material
        return None

    def get_num_of_material(self):
        return len(self.__material_list)

    def get_material_list(self):
        return self.__material_list

    def get_pattern_list(self):
        return self.__pattern_list

    def add_rectangle_pattern(self, material, position, angle, widths):
        self.__material_list.append(material)
        epsilon_type, mu_type = material.material_type
        if (epsilon_type == "tensor") or (mu_type == "tensor"):
            self.__has_tensor = True
        pattern = RectanglePatten(
            position=position,
            angle=angle,
            widths=widths,
            area= geometry.get_rectangle_area(args2)
        )
        self.__pattern_list.append(pattern)

    def add_circle_pattern(self, material, position, radius):
        self.__material_list.append(material)
        epsilon_type, mu_type = material.material_type
        if (epsilon_type == "tensor") or (mu_type == "tensor"):
            self.__has_tensor = True
        pattern = CirclePattern(
            position=position,
            radius=radius,
            area = geometry.get_circle_area(radius)
        )
        self.__pattern_list.append(pattern)

    def add_ellipse_pattern(self, material, position, angle, halfwidths):
        self.__material_list.append(Material)
        epsilon_type, mu_type = material.material_type
        if (epsilon_type == "tensor") or (mu_type == "tensor"):
            self.__has_tensor = True
        pattern = EllipsePattern(
            position=position,
            angle=angle,
            halfwidths=halfwidths
            area=geometry.get_ellipse_area(args2)
        )
        self.__pattern_list.append(pattern)

    def add_polygon_pattern(self, material, position, angle, edgepoints):
        self.__material_list.append(material)
        epsilon_type, mu_type = material.material_type
        if (epsilon_type == "tensor") or (mu_type == "tensor"):
            self.__has_tensor = True

        edge_list = np.array([point for point in edgepoints])

        pattern = PolygonPattern(
            position=position,
            angle=angle,
            edge_list=edge_list,
            area=geometry.get_polygon_area(edge_list)
        )
        self.__pattern_list.append(pattern)

    def add_grating_pattern(self, material, center, width):
        
        self.__material_list.append(material)
        epsilon_type, mu_type = material.material_type
        if (epsilon_type == "tensor") or (mu_type == "tensor"):
            self.__has_tensor = True
        pattern = GratingPattern(
            center=center,
            width=width,
            area=geometry.get_grating_area(width)
        )
        self.__pattern_list.append(pattern)

    def __is_contain_in_geometry(self, pattern1, pattern2):
        # define the center of the geometrical pattern
        position2 = np.zeros(2)
        pattern2_type = pattern2.patter_type
        if pattern2_type == "grating":
            # Quasi 1D for grating
            position2 = np.array([pattern2.center, 0])
        else:
            position2 = np.array(pattern2.position)

        # define whether the pattern is contained respectively
        pattern1_type = pattern1.pattern_type
        if pattern1_type == "grating":
            center = pattern1.center
            width = pattern1.width
            return geometry.is_contained_in_grating(center, position2[0], width)
        elif pattern1_type == "rectangle":
            position = pattern1.position
            widths = pattern1.widths
            return geometry.is_contained_in_rectangle(position, position2, widths)
        elif pattern1_type == "circle":
            position = pattern1.position
            radius = pattern1.radius
            return geometry.is_contained_in_circle(position, position2, radius)
        elif pattern1_type == "ellipse":
            position = pattern1.position
            halfwidths = pattern1.halfwidths
            return geometry.is_contained_in_ellipse(position, position2, halfwidths[0], halfwidths[1])
        elif pattern1_type == "polygon":
            position = pattern1.position
            edge_list = pattern1.edge_list
            return geometry.is_contained_in_polygon(position, position2, edge_list)

        return False

    def get_geometry_containment_relation(self):
        '''
        Obtain the containment relation between geometrical patterns
        '''
        if self.verbose:
            print("Starting to get geometry containment relation for the layer [{name}]".format(name=self.__name))

        area_list = np.array([pattern.area for pattern in self.__pattern_list])
        indexs = np.argsort(area_list)

        # Refresh the parent for each pattern based on the geometrical containment relations
        for i in range(len(area_list)):
            self.__pattern_list[indexs[i]].parent = -1
            for j in range(i,len(area_vectors)):
                pattern2 = self.__pattern_list[indexs[i]]
                pattern1 = self.__pattern_list[indexs[j]]
                if self.__is_contain_in_geometry(pattern1, pattern2):
                    self.__pattern_list[indexs[i]].parent = indexs[j]
                    break
        
        if self.verbose:
            print("End of the computations of the geometry containment relation for the layer [{name}]".format(
                name=self.__name
            ))

    def get_epsilon_matrix(self, Gx, Gy, omega_index):
        '''
        Get the Fourier transformation of the optical parameters of the layer
        '''
        nG = len(Gx)
        epsxx = np.zeros((nG, nG), dtype=complex)
        epsxy = np.zeros((nG, nG), dtype=complex)
        epsxz = np.zeros((nG, nG), dtype=complex)
        epsyx = np.zeros((nG, nG), dtype=complex)
        epsyy = np.zeros((nG, nG), dtype=complex)
        epsyz = np.zeros((nG, nG), dtype=complex)
        epszx = np.zeros((nG, nG), dtype=complex)
        epszy = np.zeros((nG, nG), dtype=complex)
        epszz = np.zeros((nG, nG), dtype=complex)

        muxx = np.zeros((nG, nG), dtype=complex)
        muxy = np.zeros((nG, nG), dtype=complex)
        muxz = np.zeros((nG, nG), dtype=complex)
        muyx = np.zeros((nG, nG), dtype=complex)
        muyy = np.zeros((nG, nG), dtype=complex)
        muyz = np.zeros((nG, nG), dtype=complex)
        muzx = np.zeros((nG, nG), dtype=complex)
        muzy = np.zeros((nG, nG), dtype=complex)
        muzz = np.zeros((nG, nG), dtype=complex)

        for i in range(len(self.__pattern_list)):
            pattern = self.__pattern_list[i]
            material = self.__material_list[i]

            # epsilon = material.get_epsilon_at_index(omega_index)
            # mu = material.get_mu_at_index(omega_index)
            material_parent = self.__background_material
            if pattern.parent != -1:
                material_parent = self.__material_list[pattern.parent]

            pattern_type = pattern.get_pattern_type()
            if pattern_type == "grating":
                dvals = self.__transform_grating(pattern, material, material_parent, Gx, Gy, omega_index)
            elif pattern_type == "rectangle":
                dvals = self.__transform_rectangle(pattern, material, material_parent, Gx, Gy, omega_index)
            elif pattern_type == "circle":
                dvals = self.__transform_circle(pattern, material, material_parent, Gx, Gy, omega_index)
            elif pattern_type == "ellipse":
                dvals = self.__transform_ellipse(pattern, material, material_parent, Gx, Gy, omega_index) 
            elif pattern_type == "polygon":
                dvals = self.__transform_polygon(pattern, material, material_parent, Gx, Gy, omega_index)

            epsxx += dvals[0]
            epsxy += dvals[1]
            epsxz += dvals[2]
            epsyx += dvals[3]
            epsyy += dvals[4]
            espyz += dvals[5]
            epszx += dvals[6]
            epszy += dvals[7]
            epszz += dvals[8]

            muxx += dvals[9]
            muxy += dvals[10]
            muxz += dvals[11]
            muyx += dvals[12]
            muyy += dvals[13]
            muyz += dvals[14]
            muzx += dvals[15]
            muzy += dvals[16]
            muzz += dvals[17]
                
        # consider the background material
        eps_background = self.__background_material.get_epsilon_at_index(omega_index)
        mu_background = self.__background_material.get_mu_at_index(omega_index)
        
        I_mat = np.identity((nG, nG), dtype=complex)
        epsxx += eps_background[0,0]*I_mat
        epsxy += eps_background[0,1]*I_mat
        epsxz += eps_background[0,2]*I_mat
        epsyx += eps_background[1,0]*I_mat
        epsyy += eps_background[1,1]*I_mat
        espyz += eps_background[1,2]*I_mat
        epszx += eps_background[2,0]*I_mat
        epszy += eps_background[2,1]*I_mat
        epszz += eps_background[2,2]*I_mat

        muxx += mu_background[0,0]*I_mat
        muxy += mu_background[0,1]*I_mat
        muxz += mu_background[0,2]*I_mat
        muyx += mu_background[1,0]*I_mat
        muyy += mu_background[1,1]*I_mat
        muyz += mu_background[1,2]*I_mat
        muzx += mu_background[2,0]*I_mat
        muzy += mu_background[2,1]*I_mat
        muzz += mu_background[2,2]*I_mat

        return [epsxx, epsxy, epsxz, 
                epsyx, epsyy, epsyz,
                epszx, epszy,  epszz,
                muxx, muxy, muxz, 
                muyx, muyy, muyz,
                muzx, muzy, muzz]
            
    def __transform_grating_element(G_mat, width):
        geometry_mat = width * geometry.sinc(G_mat / 2 * width)
        return geometry_mat

    def __transform_grating(pattern, material, material_parent, Gx, Gy, omega_index):
        center = pattern.center * MICRON
        width = pattern.width *MICRON

        Gx_r, Gx_l = np.meshgrid(Gx, Gx)
        Gmat = Gx_l - Gx_r 
        phase = np.exp(1j*Gmat*center)

        eps_parent = material_parent.get_epsilon_at_index(omega_index)
        mu_parent = material_parent.get_mu_at_index(omega_index)
        eps = material.get_epsilon_at_index(omega_index)
        mu = material.get_mu_at_index(omega_index)

        geometry_mat = self.__transform_grating_element(G_mat, width)     
        dimension = self.__lattice.get_dimension()
        if dimension == "one":
            area = self.__lattice.get_area() * MICRON
        else:
            area = self.__lattice.get_area() * np.sqaure(MICRON)

        dval = list()
        dval_mu = list()
        for i in range(3):
            for j in range(3):
                eps_mat = (eps[i,j] - eps_parent[i,j])*phase*geometry_mat
                dval.append(eps_mat)

                mu_mat = (mu[i,j] - mu_parent[i,j])*phase*geometry_mat
                dval_mu.append(mu_mat)
        
        dval.extend(dval_mu)
        return dval

    def __transform_rectangle_element(Gx_mat, Gy_mat, widths):
        widthx, widthy = widths
        geometry_mat = widthx * widthy * geometry.sinc(Gx_mat * widthx / 2) * geometry.sinc(Gy_mat * widthy /2)
        return geometry_mat

    def __transform_rectangle(pattern, material, material_parent, Gx, Gy, omega_index):
        position = pattern.position * MICRON
        angle = pattern.angle * np.pi/180
        widths = pattern.widths * MICRON

        Gx_r, Gx_l = np.meshgrid(Gx, Gx)
        Gy_r, Gy_l = np.meshgrid(Gy, Gy)
        Gx_mat = Gx_l - Gx_r
        Gy_mat = Gy_l - Gy_r
        phase = np.exp(1j*(Gx_mat*position[0] + Gy_mat*position[1]))
        G_temp = Gx_mat * np.cos(angle) + Gy_mat * np.sin(angle)
        Gy_mat = -Gx_mat * np.sin(angle) + Gy_mat * np.cos(angle)
        Gx_mat = G_temp
        geometry_mat = self.__transform_rectangle_element(Gx_mat, Gy_mat, widths)     

        eps_parent = material_parent.get_epsilon_at_index(omega_index)
        mu_parent = material_parent.get_mu_at_index(omega_index)
        eps = material.get_epsilon_at_index(omega_index)
        mu = material.get_mu_at_index(omega_index)

        dimension = self.__lattice.get_dimension()
        if dimension == "one":
            area = self.__lattice.get_area() * MICRON
        else:
            area = self.__lattice.get_area() * np.sqaure(MICRON)

        dval = list()
        dval_mu = list()
        for i in range(3):
            for j in range(3):
                eps_mat = (eps[i,j] - eps_parent[i,j])*phase*geometry_mat
                dval.append(eps_mat)

                mu_mat = (mu[i,j] - mu_parent[i,j])*phase*geometry_mat
                dval_mu.append(mu_mat)
        
        dval.extend(dval_mu)
        return dval

    def __transform_circle_element(Gx_mat, Gy_mat, radius):
        rho = np.sqrt(np.sqaure(Gx_mat) + np.sqaure(Gy_mat)) * radius
        jinc_mat = geometry.jinc(rho)
        geometry_mat = 2*np.pi*np.square(radius) * jinc_mat
        return geometry_mat

    def __transform_circle(pattern, material, material_parent, Gx, Gy, omega_index):
        position = pattern.position * MICRON
        radius = pattern.radius * MICRON

        Gx_r, Gx_l = np.meshgrid(Gx, Gx)
        Gy_r, Gy_l = np.meshgrid(Gy, Gy)
        Gx_mat = Gx_l - Gx_r
        Gy_mat = Gy_l - Gy_r
        phase = np.exp(1j*(Gx_mat*position[0] + Gy_mat*position[1]))
        geometry_mat = self.__transform_circle_element(Gx_mat, Gy_mat, radius)     

        eps_parent = material_parent.get_epsilon_at_index(omega_index)
        mu_parent = material_parent.get_mu_at_index(omega_index)
        eps = material.get_epsilon_at_index(omega_index)
        mu = material.get_mu_at_index(omega_index)

        dimension = self.__lattice.get_dimension()
        if dimension == "one":
            area = self.__lattice.get_area() * MICRON
        else:
            area = self.__lattice.get_area() * np.sqaure(MICRON)

        dval = list()
        dval_mu = list()
        for i in range(3):
            for j in range(3):
                eps_mat = (eps[i,j] - eps_parent[i,j])*phase*geometry_mat
                dval.append(eps_mat)

                mu_mat = (mu[i,j] - mu_parent[i,j])*phase*geometry_mat
                dval_mu.append(mu_mat)
        
        dval.extend(dval_mu)
        return dval

    def __transform_ellipse_element(Gx_mat, Gy_mat, halfwidths):
        a, b = halfwidths
        rho = np.square(np.sqaure(a*Gx_mat) + np.square(b*Gy_mat))
        jinc_mat = geometry.jinc(rho)
        geometry_mat = 2*np.pi*a*b*jinc_mat
        return geometry_mat

    def __transform_ellipse(pattern, material, material_parent, Gx, Gy, omega_index):
        position = pattern.position * MICRON
        angle = pattern.angle * np.pi/180
        halfwidths = pattern.halfwidths * MICRON

        Gx_r, Gx_l = np.meshgrid(Gx, Gx)
        Gy_r, Gy_l = np.meshgrid(Gy, Gy)
        Gx_mat = Gx_l - Gx_r
        Gy_mat = Gy_l - Gy_r
        phase = np.exp(1j*(Gx_mat*position[0] + Gy_mat*position[1]))
        G_temp = Gx_mat * np.cos(angle) + Gy_mat * np.sin(angle)
        Gy_mat = -Gx_mat * np.sin(angle) + Gy_mat * np.cos(angle)
        Gx_mat = G_temp
        geometry_mat = self.__transform_ellipse_element(Gx_mat, Gy_mat, halfwidths)     

        eps_parent = material_parent.get_epsilon_at_index(omega_index)
        mu_parent = material_parent.get_mu_at_index(omega_index)
        eps = material.get_epsilon_at_index(omega_index)
        mu = material.get_mu_at_index(omega_index)

        dimension = self.__lattice.get_dimension()
        if dimension == "one":
            area = self.__lattice.get_area() * MICRON
        else:
            area = self.__lattice.get_area() * np.sqaure(MICRON)

        dval = list()
        dval_mu = list()
        for i in range(3):
            for j in range(3):
                eps_mat = (eps[i,j] - eps_parent[i,j])*phase*geometry_mat
                dval.append(eps_mat)

                mu_mat = (mu[i,j] - mu_parent[i,j])*phase*geometry_mat
                dval_mu.append(mu_mat)
        
        dval.extend(dval_mu)
        return dval

    def __transform_polygon_element(Gx_mat, Gy_mat, edge_list):
        polygon_area = geometry.get_polygon_area(edge_list)
        geometry_mat = np.zeros(Gx_mat.shape, dtype=complex)

        nGx, nGy = Gx_mat.shape
        for i in range(nGx):
            for j in range(nGy):
                u = Gx_mat[i,j]
                v = Gy_mat[i,j]

                if (u == 0) and (v == 0):
                    geometry_mat[i,j] = polygon_area
                else:
                    length_edge_list = len(edge_list)
                    for index in range(length_edge_list):
                        x_cur, y_cur = edge_list[index]
                        if (index == length_edge_list - 1):
                            x_next, y_next = edge_list[0]
                        else:
                            x_next, y_next = edge_list[i+1]

                if (u == 0):
                    val = 1j/v*(x_next - x_cur) * np.exp(1j*(u*(x_next + x_cur)/2 + v*(y_next + y_cur)/2)) * geometry.sinc((x_next - x_cur)*u /2 + (y_next-y_cur)*v/2)
                else:
                    val = -1j/u*(y_next - y_cur) * np.exp(1j*(u*(x_next + x_cur)/2 + v*(y_next + y_cur)/2)) * geometry.sinc((x_next - x_cur)*u /2 + (y_next-y_cur)*v/2)
                
                geometry_mat[i,j] += val

        return geometry_mat

    def __transform_polygon(pattern, material, material_parent, Gx, Gy, omega_index):
        position = pattern.position * MICRON
        angle = pattern.angle * np.pi/180
        edge_list = pattern.edge_list * MICRON

        Gx_r, Gx_l = np.meshgrid(Gx, Gx)
        Gy_r, Gy_l = np.meshgrid(Gy, Gy)
        Gx_mat = Gx_l - Gx_r
        Gy_mat = Gy_l - Gy_r
        phase = np.exp(1j*(Gx_mat*position[0] + Gy_mat*position[1]))
        G_temp = Gx_mat * np.cos(angle) + Gy_mat * np.sin(angle)
        Gy_mat = -Gx_mat * np.sin(angle) + Gy_mat * np.cos(angle)
        Gx_mat = G_temp
        geometry_mat = self.__transform_polygon_element(Gx_mat, Gy_mat, edge_list)     

        eps_parent = material_parent.get_epsilon_at_index(omega_index)
        mu_parent = material_parent.get_mu_at_index(omega_index)
        eps = material.get_epsilon_at_index(omega_index)
        mu = material.get_mu_at_index(omega_index)

        dimension = self.__lattice.get_dimension()
        if dimension == "one":
            area = self.__lattice.get_area() * MICRON
        else:
            area = self.__lattice.get_area() * np.sqaure(MICRON)

        dval = list()
        dval_mu = list()
        for i in range(3):
            for j in range(3):
                eps_mat = (eps[i,j] - eps_parent[i,j])*phase*geometry_mat
                dval.append(eps_mat)

                mu_mat = (mu[i,j] - mu_parent[i,j])*phase*geometry_mat
                dval_mu.append(mu_mat)
        
        dval.extend(dval_mu)
        return dval

class Structure:
    '''
    Define the whole multilayer structure composed of layers
    '''
    def __init__(self, lattice, material_list, layer_list, verbose=False):
        self.__layer_list = layer_list
        self.__material_list = material_list
        self.__lattice = lattice

    def __init__(self, structure):
        layer_list = structure.get_layer_list()
        material_list = structure.get_material_list()
        lattice = structure.get_lattice()
        self.__layer_list = layer_list.copy()
        self.__material_list = material_list.copy()
        self.__lattice = lattice

    def get_lattice(self):
        return self.__lattice

    def add_material(self, material):
        self.__material_list.append(material)

    def get_material_list(self):
        return self.__material_list

    def add_layer(self, layer):
        self.__layer_list.append(layer)

    def set_lattice(self, lattice):
        self.__lattice = lattice

    def delete_layer_by_name(self, name):
        index = None
        for i in range(len(self.__layer_list)):
            layer = self.__layer_list[i]
            if layer.name == name:
                index = i
                break
        if index is not None:
            del self.__layer_list[index]

    def delete_layer_by_layer(self, layer):
        self.delete_layer_by_name(layer.name)

    def get_layer_by_index(self, index):
        if index < 0:
            raise Exception("Index: smaller than zero for the layer list")
        if len(self.__layer_list) <= index:
            #return None
            raise Exception("Index: larger than the length of the layer list")
        return self.__layer_list[index]    

    def get_layer_by_name(self, name):
        for layer in self.__layer_list:
            if layer.name == name:
                return layer
        return None

    def get_num_of_layer(self):
        return len(self.__layer_list)

    def get_thickness_list(self):
        thickness_list = list()
        for layer in self.__layer_list:
            thickness_list.append(layer.get_thickness())
        return thickness_list

    def get_layer_list(self):
        return self.__layer_list.copy()

    def delete_layer_by_index(self, index):
        if index < 0:
            raise Exception("Index: smaller than zero for the layer map")
        if len(self.__layer_list) <= index:
            raise Exception("Index: larger than the length of the layer map")
        del self.__layer_list[index]   


    