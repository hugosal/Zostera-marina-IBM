#o-*- coding: utf-8 -*-
#$Id: $

#functions to create  a sample world for the model

import os
import math
import pickle
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np

def world(n, m):
    #this function is to create the grid of x and y coordinates of a world
    #the grid will be of evenly space coordinates from -n/2 to n/2, so n  and m are the total
    #widht and height of the world, the world is made of 50 *50 cells, this can be changed adjusting
    #in the following steps.
    #to create the grid use the meshgrid, where the number of interval between the max and minumun a determined by 
    #the complex number j, in this case 50
    grid_x, grid_y = np.mgrid[-n/2:n/2:50j, -m/2:m/2:50j] 
    #in this case, also create a  batimetry of the world by appying a depth fuction to the information on the
    depth_x_y = depth(grid_x, grid_y) # the depth fuction is defined elsewhere
    return (grid_x, grid_y, depth_x_y) #return all three components of the world

#the depth function to create a world shaped like a parabolic canal, with a bottom of -4 at the
#middleof th canalthe depth function is
def depth(x,y):
    return (0.0000009*(y ** 2)-4000)

#for a circular 'hole' world the depth fuction would be
#def depth(x, y):
#    return (x ** 2 + y ** 2) / 800

#to create and save a world to use in the simulation:
grid_x, grid_y,depth_x_y = world(200000, 200000)
place = open(os.getcwd() + os.sep + 'cannal_200m_broad_4m_prof.dat', "wb")
pickle.dump((grid_x, grid_y,depth_x_y ), place)
place.close()
#now you can use the canal world in a simulation
print("Succes")
