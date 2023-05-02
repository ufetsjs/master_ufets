# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 14:01:24 2021

@author: xy0264
"""
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm

plt.close('all')

# General parameters
NX=20
NY=20
VX=1
VY=1
DT=0.1
T_MAX=100
NUM_T=int(round(T_MAX/DT))
C_IN=np.zeros(NUM_T+1)
C_IN[0]=1
C_IN=np.ones(NUM_T+1)*1

# Export frames?
exp_frame=True
T_SKIP=round(NUM_T/20)

# Minimum for plot
C_MIN=VX*np.min(C_IN)*DT*1e-2
    
# Setup Concentration matrix
c=np.zeros((NX,NY,NUM_T+1))

# Time loop
for t in range(NUM_T):
       
    # Initialize concentration change
    dc=np.zeros((NX,NY))
    
    # Go through all compartments
    for x in range(NX):
        for y in range(NY):
            # In the inner matrix
            if (x != 0) and (y != 0) and (y != NY-1):
                dc[x,y]=c[x-1,y,t]*VX+c[x,y-1,t]*VY-c[x,y,t]*VX-c[x,y,t]*VY
            # RB links
            if (x == 0) and (y != 0) and (y != NY-1):
                dc[x,y]=C_IN[t]*VX+c[x,y-1,t]*VY-c[x,y,t]*VX-c[x,y,t]*VY
            if (x != 0) and (y == NY-1):
                dc[x,y]=c[x-1,y,t]*VX+c[x,y-1,t]*VY-c[x,y,t]*VX
            if (x != 0) and (y == 0):
                dc[x,y]=c[x-1,y,t]*VX-c[x,y,t]*VX-c[x,y,t]*VY
            if x == 0 and y == 0:
                dc[x,y]=C_IN[t]*VX-c[x,y,t]*VX-c[x,y,t]*VY
            if x == 0 and y == NY-1:
                dc[x,y]=C_IN[t]*VX-c[x,y,t]*VX+c[x,y-1,t]*VY
    
    #print(dc)
    # Update concentration matrix
    c[:,:,t+1]=c[:,:,t]+dc*DT

# Plot setup
fig = plt.figure()
ax = fig.add_subplot()

# Set zeros in c to C_MIN
c[c==0]=C_MIN

# Time loop
for t in range(0,NUM_T,T_SKIP):
    # Clear axes
    ax.cla()
    
    cp = ax.pcolormesh(np.arange(NX),np.arange(NY),c[:,:,t].T,norm=LogNorm(vmin=C_MIN, vmax=np.max(c)),edgecolor='None',shading='nearest')
    # print('cp' in locals())
    # # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    # if 'cb' in locals():
    #     cp.remove()
    cp = plt.colorbar(cp,cax)
    cp.set_label('particle concentration $c$')
    
    # if exp_frame:
    #     plt.savefig(f"frames\gif_{t:05d}.png",dpi=200)
        
    # plt.pause(DT)

plt.show()
    


