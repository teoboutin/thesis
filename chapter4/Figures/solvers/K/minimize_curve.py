







import random
import math
import numpy as np
import matplotlib.pyplot as plt

import copy
def dist(p1, p2):
    return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2
    
    
# doesn't work. actually might if I find a better way to insert on correct side
# but the problem is much simpler to solve (and much harder to, because it is a NP hard problem, I am very rusted...)
def minimize_curve_length(list_of_points):
    # ~ print(list_of_points)
    
    result=[list_of_points[0]]
    for point_to_insert in list_of_points[1:]:
        # ~ print(result)
        # find the position of closest point:
        min_val=np.infty
        min_pos=0
        for i in range(len(result)):
            pp=result[i]
            d=dist(point_to_insert,pp)
            if d<min_val:
                min_val=d
                min_pos=i
        
        # insert the point on the correct side:
        closest=result[min_pos]
        if len(result)==1:
            result.insert(0, point_to_insert)
        elif min_pos==0:
            if dist(point_to_insert, result[1]) > dist(result[0], result[1]): 
                # point we want to insert is farther to second point than first point: then insert outside
                result.insert(0, point_to_insert)
            else:
                # point is between the two first points, insert between
                result.insert(1, point_to_insert)
        elif min_pos == len(result)-1:
            n=len(result)
            if dist(point_to_insert, result[n-2]) > dist(result[n-1], result[n-2]): 
                # point we want to insert is farther to second last point than last point: then insert outside
                result.insert(min_pos+1, point_to_insert)
            else:
                # point is between the two first points, insert between
                result.insert(min_pos, point_to_insert)
        else:
            dist_left=dist(point_to_insert, result[min_pos-1])
            dist_right=dist(point_to_insert, result[min_pos+1])
            if dist_left < dist_right:
                result.insert(min_pos, point_to_insert)
            else:
                result.insert(min_pos+1, point_to_insert)
    return result

# this solution is probably optimal if on the contour of a convex shape... kinda restrictive. hyperbolas won't work well

def closest(ref, lop):
    min_val=np.infty
    min_pos=0
    for i in range(len(lop)):
        pp=lop[i]
        d=dist(ref,pp)
        if d<min_val:
            min_val=d
            min_pos=i
    return min_pos, min_val
    
def solve_TSP(list_of_points, start=0):
    lop=copy.deepcopy(list_of_points)
    
    result=[lop.pop(start)]
    total_dist=0
    
    while len(lop)>0:
        
        # search closest point to each extremety. 
        
        min_pos_left, min_val_left  = closest(result[0], lop)
        min_pos_right, min_val_right  = closest(result[-1], lop)
        
        # then append the closest one
        # better handling of non closed curves
        # still bad for disjoint spaces (obviously)
        if min_val_left < min_val_right:
            c = lop.pop(min_pos_left)
            total_dist+=min_val_left
            result.insert(0,c)
        else:
            c = lop.pop(min_pos_right)
            total_dist+=min_val_right
            result.append(c)
    return result, total_dist
    
def compute_length(lop):
    total_dist=0
    for i in range(1,len(lop)):
        total_dist+=dist(lop[i-1], lop[i])
    return total_dist
    
def solve_TSP_not_closed(list_of_points):
    result, best_dist=solve_TSP(list_of_points, 0)
    best_result= copy.deepcopy(result)
    # shift the values
    for i in range(0,len(list_of_points)):
        last=result.pop()
        result.insert(0,last)
        d = compute_length(result)
        if d < best_dist:
            best_dist=d
            best_result = copy.deepcopy(result)
        
            
    return best_result


def plot_list_of_pts(lp, **kwargs):
    x = [p[0] for p in lp]
    y = [p[1] for p in lp]
    plt.plot(x,y, **kwargs)

if __name__ == "__main__":
    
    n=100
    
    x = np.random.rand(n)

    def fun(x):
        return (x-0.5)**2 + 1
    y = fun(x)
    
    
    testlist=[(x[i],y[i]) for i in range(len(x))]

    
    plot_list_of_pts(testlist, label="test", marker="o")
    ts = solve_TSP_not_closed(testlist)
    plot_list_of_pts(ts, label="sorted", marker="x")
    
    plt.legend()
    plt.show()