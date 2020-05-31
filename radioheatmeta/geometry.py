import numpy as np 

def get_rectangle_area(args):
    return args[0]*args[1]

def get_ellipse_area(args):
    return np.pi*args[0]*args[1]

def get_circle_area(radius):
    return np.pi*radius*radius

def cross_product(point1, point2):
    x1, y1 = point1
    x2, y2 = point2
    return x1*y2 - y1*x2

def get_polygon_area(edge_list):
    area = 0
    for i in range(len(edge_list)-1):
        area += cross_product(edge_list[i], edge_list[i+1])
    area += cross_product(edge_list[-1], edge_list[0])
    return np.abs(area)/2

def get_grating_area(x):
    return x