#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 14:54:45 2020

@author: goharshoukat
"""
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.tri as tri
from pydmd import DMD
import scipy
from mpl_toolkits.mplot3d.axes3d import Axes3D
from decimal import Decimal
import math


cmMotion_data = np.loadtxt('cmMotion')
#drop zero columns from the cmMotion
cmMotion_data = np.delete(cmMotion_data, [1,2],1)
#data slicing to extract only relevant time steps from  30s <= t <= 50s

#following operations are to combine the p and v data with CoM. 

#the CM file contains redundant time steps. the data writen for P and v isnt 
#in those particular denominations and hence further data extraction has to be
#carried out. 7th iteration is the data we need to keep. the rest can be ignored. 

cmMotion = ([])
for i in range(np.where(cmMotion_data[:,0] == 30)[0][0],len(cmMotion_data)):
    if Decimal(str(cmMotion_data[i,0]))%Decimal('0.01') == Decimal('0.00'):
       cmMotion = np.append(cmMotion,cmMotion_data[i,1])
#transpose cmMotion to make it compatible with P_v_data file. 
cmMotion = cmMotion.reshape(-1,len(cmMotion))
#read all timesteps with information about velocity, and pressure. 
#the file contains velocity vectors and pressure
path = os.path.join('DMD')
#to read in a sequence
#all file path names transferred. 
all_files = sorted(glob.glob(os.path.join(path, '*'))) 

#pressure and velocity data read
p_v_data =[]
for file in all_files:
    data_array = np.loadtxt(file)
    p_v_data = np.append(p_v_data,data_array)
#resize p_v_data to divide into 128 rows and 5000 time steps. 
p_v_data = p_v_data.reshape(len(all_files),-1).T

len(p_v_data)
#Append the two arrays p_v_data and cmMotion
p_v_cm = np.append(p_v_data,cmMotion,0)
#combine the data for CoM with the pressure and velocity

##############################################################################
#creating a useful x-axis measurement to plot the data against for visualsization
##############################################################################

r = 0.2 #radius of cylinder around which measurements are taken
#radius is not the surface of the structure due to no-slip condition. 
#offset applied

#number of points along which measurements are taken
vertices = 32

circumference = 2 * math.pi * r
#open up the circle, assume it to be straight line
# and divide circumference into 32 points. 
circ_delta = circumference/vertices
#create an array
circ_x = np.arange(0, circumference,circ_delta)

#time array
time_start = 30
time_end = 50
time_delta = 0.01
time = np.arange(time_start, time_end+time_delta, time_delta)
time.T
xs, ts = np.meshgrid(circ_x, time)
#actual data
#following 3 arrays are redundant but increase readability
act_press = p_v_data[0::4]
act_vel_x = p_v_data[1::4]
act_vel_z = p_v_data[3::4]
act_CM = p_v_cm[-1]

#plt.plot(circ_x,act_press[:,0])
'''
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(xs, ts, act_press.T, rstride=1, cstride=1, cmap='hot')
plt.show()

fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(xs, ts, act_vel_x.T, rstride=1, cstride=1, cmap='hot')
ax.plot_surface(xs, ts, reconst_vel_x.T, rstride=1, cstride=1, cmap='cold')
plt.show()
'''
##############################################################################
                                    #DMD begins
##############################################################################
dmd = DMD()
decomposed_data = dmd.fit(p_v_cm)
dmd.plot_eigs()

modes = dmd.modes
vel_x_modes = modes[1::4]
#plot few modes to represent shape of modes. It becomes too clustered otherwise. 
plt.plot(circ_x, vel_x_modes[:,0:3])    



for dynamic in dmd.dynamics:
    plt.plot(dynamic.real)
    plt.title('Dynamics')
plt.show()

reconstructed_data = dmd.reconstructed_data
##following 3 arrays are redundant but increase readability
reconst_press = reconstructed_data[0:len(p_v_data):4]
reconst_vel_x = reconstructed_data[1::4]
reconst_vel_z = reconstructed_data[3::4]
reconst_CM = reconstructed_data[-1]


##############################################################################
                            #####write images to files#####
##############################################################################
if not os.path.exists('vel_x'):
    os.makedirs('vel_x')

if not os.path.exists('vel_z'):
    os.makedirs('vel_z')

if not os.path.exists('pressure'):
    os.makedirs('pressure')
    

if not os.path.exists('CM'):
    os.makedirs('CM')
'''
fig = plt.figure()
plt.plot(time, act_CM, label= 'Actual')
plt.plot(time, reconst_CM, label = 'Reconstructued')
plt.title('Center of Mass motion')
plt.legend()
plt.xlabel('Time [s]')
plt.ylabel('Z coordinate[m]')
fig.savefig('CM.png')
plt.close(fig)

#save files in a folder for comparison with velocity x component
path_vel_x = os.path.join('vel_x')
for i in range(len(all_files)):
    new_file = str(i) + '.png'
    fig1 = plt.figure()
    plt.plot(circ_x,act_vel_x[:,i], label='Actual Data')
    plt.plot(circ_x, reconst_vel_x[:,i], label='Reconstructed Data')
    plt.title('X compment of Velocity (m/s)')
    plt.xlabel('Circumference(m)')
    plt.ylabel('Velocity (m/s)')
    plt.legend()
    fig1.savefig(os.path.join(path_vel_x,new_file))
    plt.close(fig1)


#save files in a folder for comparison with velocity z component
path_vel_z = os.path.join('vel_z')
for i in range(len(all_files)):
    new_file = str(i) + '.png'
    fig2 = plt.figure()
    plt.plot(circ_x,act_vel_z[:,i], label='Actual Data')
    plt.plot(circ_x, reconst_vel_z[:,i], label='Reconstructed Data')
    plt.title('Z compment of Velocity (m/s)')
    plt.xlabel('Circumference(m)')
    plt.ylabel('Velocity (m/s)')
    plt.legend()
    fig2.savefig(os.path.join(path_vel_z,new_file))
    plt.close(fig2)
    
#save files in a folder for comparison with velocity z component
path_press = os.path.join('pressure')
for i in range(len(all_files)):
    new_file = str(i) + '.png'
    fig3 = plt.figure()
    plt.plot(circ_x,act_press[:,i], label='Actual Data')
    plt.plot(circ_x, reconst_press[:,i], label='Reconstructed Data')
    plt.title('Pressure')
    plt.xlabel('Circumference(m)')
    plt.ylabel('Pressure[Pa]')
    plt.legend()
    fig3.savefig(os.path.join(path_press,new_file))
    plt.close(fig3)

##############################################################################
'''