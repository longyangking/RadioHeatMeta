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

def is_contained_in_grating(center1, center2, width):
    return np.abs(center2 - center1) <= width/2

def is_contained_in_rectangle(center1, center2, width):
    x1, y1 = center1
    x2, y2 = center2
    dx, dy = width
    return (np.abs(x2 - x1) <= dx/2) and (np.abs(y2 - y1) <= dy/2)

def is_contained_in_circle(center1, center2, radius):
    center1 = np.array(center1)
    center2 = np.array(center2)
    return np.sum(np.square(np.abs(center1 - center2))) <= np.square(np.abs(radius))

def is_contained_in_ellipse(center1, center2, a, b):
    center1 = np.array(center1)
    center2 = np.array(center2)
    r = np.array([a,b])
    return np.sum(np.square(np.abs(center1 - center2)/r)) <= 1

def is_contained_in_polygon(center1, center2, edge_list):
    is_inside = False
    center1 = np.array(center1)
    center2 = np.array(center2)
    r = center2 - center1

    # Determine is_contained based on vector operations
    for i in range(len(edge_list)):
        s1 = (edge_list[i][1] > r[1]) != (edge_list[i-1][1] > r[1])
        s2 = r[0] < (edge_list[i-1][0] - edge_list[i][0])*(r[1] - edge_list[i][1])/(edge_list[i-1][1] - edge_list[i][1]) + edge_list[i][0]
        if s1 and s2:
            is_inside = not is_inside

    return is_inside

def on_segment(point1, point2, porint3):
    # TODO
    pass

def orientation(point1, point2, point3):
    # TODO
    pass

def do_intersect(p1, q1, p2, q2):
    # TODO
    pass