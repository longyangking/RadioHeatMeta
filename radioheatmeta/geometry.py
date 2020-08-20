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

#  Define several geometry operations for structures

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

def on_segment(p, q, r):
    if (q[0] <= np.max([p[0], r[0]])) and \
        (q[0] >= np.min([p[0], r[0]])) and \
        (q[1] <= np.max([p[1], r[1]])) and \
        (q[1] >= np.min([p[1], r[1]])):
        return True
    return False

def orientation(p, q, r):
    val = (q[1] - p[1])*(r[0] - q[0]) - (q[0] - p[0])*(r[1] - q[1])
    if (val == 0):
        return 0
    if val > 0:
        return 1
    else:
        return 2

def do_intersect(p1, q1, p2, q2):
    o1 = orientation(p1, q1 ,p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    if (o1 != o2) and (o3 != o4):
        return True
    if (o1 == 0) and on_segment(p1, p2, q1):
        return True
    if (o2 == 0) and on_segment(p1, q2, q1):
        return True
    if (o3 == 0) and on_segment(p2, p1, q2):
        return True
    if (o4 == 0) and on_segment(p2, q1, q2):
        return True
    
    return False

def sinc(x):    # TODO Transform into the matrix operation
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
    if x == 0.0:
        return 0.5
    return sp.special.j1(x)/x