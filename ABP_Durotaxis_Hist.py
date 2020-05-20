#!/usr/bin/env python
# coding: utf-8

# ## Simulating  Active Brownian motion
# In this notebook, we simulate the trajectory of a single, non-interacting active Brownian particle. Unlike passive Brownian particles (e.g. pollen grains in water) that do diffusive random walks, active Brownian particles (ABPs) propel themselves with directed motion in addition to diffusing randomly. 
# 
# The APB model can be used to accurately desribe the motion of certain types of cells, those which locomote actively using focal adhesions. The distribution of focal adhesions within a cell may vary continuously, tending to reorient the cell with respect to it's past momentum and giving rise to different speeds with which the cell moves. It is this angular reorientation which can be modelled as a diffusive process.
# 
# Durotaxis is the preference for cells to move preferentially towards a stiffer substrate. In this project we implement durotaxis by incorporating the dependence of a cell's rotational diffusion properties on the stiffness of the medium upon which a cell is crawling.
# 
# *It should be noted that although all of the simulations in this project display the trajectories or data associated with numerous cells, intercellular interactions are not taken into account here. This is a simple single cell model.*

# In[ ]:


import argparse

parser = argparse.ArgumentParser()

parser.add_argument('L', default = 3, help = 'width of stiffness transition region', type = float)

parser.add_argument('Dr', default = 5, help = 'baseline rotational diffusion constant', type = float)

parser.add_argument('ks', default = 1, help = 'stiffness on soft side of system', type = float)

parser.add_argument('kh', default = 50, help = 'stiffness on hard side of system', type = float)

parser.add_argument('vel', default = 1, help = 'inherent cellular velocity', type = float)

parser.add_argument('nw', default = 1000, help = '# of cells used to collect statistical data', type = int)

parser.add_argument('bb', default = 20, help = 'size of system', type = int)

args = parser.parse_args()

L = args.L

Drot0 = args.Dr

ks = args.ks

kh = args.kh 

vel = args.vel

num_walks = args.nw

Bbox = args.bb


# In[3]:


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numpy.random import random as rand
from scipy import stats

import os.path
from os import path


# ### 2D Confined ABP model with durotaxis: an overview.
# Cells perform a ballistic step (proportional to $dt$) and random walk in the space of $\theta$-values. 
# 
# Furthermore we implement walled boundary conditions, treating a wall as a reflecting surface. If a walker goes beyond the confining boundary of the box, it's position is altered so as to keep it in the box. Note that when this reflection of a walker by a boundary occurs, the velocity of the walker is not infuenced by the reflection.
# 
# We also incorporate 1D substrate stiffness into the model. Specifically the stiffness will depend on the location of a walker within the confining box. To the left (-x) the substrate is more soft, to the right (+x) it is more stiff, and in the center is a transition region. The rotational diffusion constant of each walker is then made to depend on the stiffness of the substrate upon which it is walking, so that a walker in the stiffer region will be taking steps more ballistically, while a walker in the softer region will be taking steps more diffusively. 
# 
# To extract meaningful data, we collect statistics from an ensemble of identically prepared systems, i.e. individual cells which are placed in identical environments, whose ballistic and diffusive components of motion are drawn from the same distributions.

# # Statistical data alongside trajectory plots
# In the following cell we plot the trajectories of ensembles of cells within systems with varying substrate stiffness.
# 
# We also extract the histogram data corresponding to the horizontal endpoints of each cell.

# In[5]:


def truncate(num, decimal = 0):
    mult = 10 ** decimal
    return int(num * mult)/mult

#ABP model

#ABP model parameters: ballistic velocity, time step, rotational diffusion constant, translational diffusion constant
dt = 0.01; Dtrans = 0.001;

#outer frame used for plot
Bframe = Bbox/10

plt.figure(figsize=(10,5))
    
xdata = np.empty(num_walks)
    
plt.subplot(121)
    
#loop through cells in the ensemble
for j in range(num_walks):
    
    # initialize arrays that store x,y and theta values, as well as initial particle position and angle

    xvec=np.zeros(0); yvec=np.zeros(0); thetavec = np.zeros(0)
    x=0.0; y = 0.0; theta = (2*np.pi)*rand(1)
    xvec = np.append(xvec,x); yvec = np.append(yvec,y); thetavec = np.append(thetavec, theta)

    #number of discrete steps each cell takes
    num_steps = 1000;
        
    #loop through steps for each cell
    for i in range(num_steps):
        
        #determine substrate stiffness based on position
        if x > -L and x < L:
            k = ks + ((kh-ks)/(2*L))*(x + L)
        elif x > L:
            k = kh
        elif x < -L:
            k = ks
            
        #calculate rotational diffusion constant based on substrate stiffness
        Drot = Drot0/k
            
        # calculate diffusive/random steps. For the x- and y-,we generate 
        #a random number between -1 & 1 to ensure that the walker can step in both directions(up/down and left/right).
        dx = np.sqrt(2*Dtrans*dt)*2*(rand(1)-0.5); 
        dy= np.sqrt(2*Dtrans*dt)*2*(rand(1)-0.5); 
        dtheta = np.sqrt(2*Drot*dt)*(2*np.pi)*(rand(1) - 0.5);
        # update coordinates (including ballistic step)
        x += vel*dt*np.cos(theta) + dx 
        y += vel*dt*np.sin(theta) + dy
        # implement walled boundary conditions
        if x > Bbox/2:
            x -= 2*(x - Bbox/2)
        elif x < -Bbox/2:
            x -= 2*(x + Bbox/2)
        elif x < Bbox/2 and x > -Bbox/2:
            x += 0
            
        if y > Bbox/2:
            y -= 2*(y - Bbox/2)
        elif y < -Bbox/2:
            y -= 2*(y + Bbox/2)
        elif y < Bbox/2 and y > -Bbox/2:
            y += 0
            
        
        # store successive positions in arrays
        xvec = np.append(xvec,x); yvec = np.append(yvec,y) 
        # update the angle and store in array
        theta += dtheta
        thetavec = np.append(thetavec, theta)
    
    #store horizontal endpoint of the jth cell
    xdata[j] = xvec[-1]
    
    #only plot trajectories of 20 cells (1000 trajectories is too much data to visualize on one plot)
    if j < 20:
        plt.plot(xvec,yvec, '-')
        plt.plot(xvec[-1], yvec[-1], 'k.')
        plt.axis([-Bbox/2 - Bframe ,Bbox/2 + Bframe,-Bbox/2 - Bframe,Bbox/2 + Bframe], 'equal');
        plt.xlabel('x')
        plt.ylabel('y')
        plt.gca().set_aspect(1)
plt.axvline(x = L, linewidth = 0.5, color = 'k', linestyle = '--');
plt.axvline(x = -L, linewidth = 0.5, color = 'k', linestyle = '--');
    
xmean = truncate(np.mean(xdata), 3)
xstd = truncate(np.std(xdata), 3)
    
plt.subplot(122)
plt.hist(xdata,20)
plt.title(f'mean = {xmean}, stdev = {xstd}')
    
plt.savefig(f'dtax_traj+hist, L ={L}.png')


# The vertical dashed lines on the trajectory plot enclose the transition region.

# In[ ]:




