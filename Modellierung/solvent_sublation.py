# -*- coding: utf-8 -*-
"""
Created on Tue May  2 09:35:48 2023

@author: jansc
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm

class solvent_sublation:
    def __init__(self):
        self.number_x = 20
        self.number_y = 20
        self.volume_flux_x = 1
        self.volume_flux_y = 1
        self.delta_t = 0.1
        self.max_t = 100
        self.number_timesteps = int(round(self.max_t / self.delta_t))
        self.concentration_in = np.zeros(self.number_timesteps + 1)
        self.concentration_in[0] = 1
        self.concentration_in=np.ones(self.number_timesteps + 1) * 1

        # Export frames?
        exp_frame = True
        T_SKIP = round(self.number_timesteps / 20)

        # Minimum for plot
        self.concentration_minimum = self.volume_flux_x * np.min(self.concentration_in) * self.delta_t * 1e-2
            
        # Setup Concentration matrix
        self.concentration = np.zeros((self.number_x, self.number_y, self.number_timesteps + 1))
        self.run_flotation()
        self.plot_certain_timesteps(range(0, self.number_timesteps, T_SKIP), exp_frame)

    def run_flotation(self):
        for t in range(self.number_timesteps):
               
            # Initialize concentration change
            dc = np.zeros((self.number_x, self.number_y))
            
            # Go through all compartments
            for x in range(self.number_x):
                for y in range(self.number_y):
                    # In the inner matrix
                    if (x != 0) and (y != 0) and (y != self.number_y - 1):
                        dc[x, y] = self.concentration[x - 1, y, t] * self.volume_flux_x + self.concentration[x, y-1 ,t] * self.volume_flux_y - self.concentration[x, y, t]*self.volume_flux_x - self.concentration[x, y, t]*self.volume_flux_y
                    # RB links
                    if (x == 0) and (y != 0) and (y != self.number_y - 1):
                        dc[x, y]=self.concentration_in[t]*self.volume_flux_x + self.concentration[x, y - 1, t] * self.volume_flux_y - self.concentration[x, y, t] * self.volume_flux_x - self.concentration[x, y, t] * self.volume_flux_y
                    if (x != 0) and (y == self.number_y - 1):
                        dc[x, y] = self.concentration[x - 1, y, t] * self.volume_flux_x + self.concentration[x, y-1, t] * self.volume_flux_y - self.concentration[x, y, t] * self.volume_flux_x
                    if (x != 0) and (y == 0):
                        dc[x, y] = self.concentration[x - 1, y, t] * self.volume_flux_x - self.concentration[x, y, t] * self.volume_flux_x - self.concentration[x, y, t] * self.volume_flux_y
                    if x == 0 and y == 0:
                        dc[x, y] = self.concentration_in[t] * self.volume_flux_x - self.concentration[x, y, t] * self.volume_flux_x - self.concentration[x, y, t] * self.volume_flux_y
                    if x == 0 and y == self.number_y - 1:
                        dc[x, y] = self.concentration_in[t] * self.volume_flux_x - self.concentration[x, y, t] * self.volume_flux_x + self.concentration[x, y - 1, t] * self.volume_flux_y
            
            #print(dc)
            # Update concentration matrix
            self.concentration[:, :, t+1] = self.concentration[:, :, t] + dc * self.delta_t
            

    def plot_certain_timesteps(self, timesteps, save):
        # Plot setup
        fig = plt.figure()
        ax = fig.add_subplot()

        # Set zeros in c to C_MIN
        self.concentration[self.concentration == 0] = self.concentration_minimum

        # Time loop
        for t in timesteps:
            # Clear axes
            ax.cla()
            cp = ax.pcolormesh(np.arange(self.number_x), np.arange(self.number_y),
                               self.concentration[:, :, t].T,norm = LogNorm(vmin = self.concentration_minimum, vmax = np.max(self.concentration)),
                               edgecolor = 'None', shading = 'nearest')
            # Colorbar
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            cp = plt.colorbar(cp, cax)
            cp.set_label('particle concentration $c$')
            # Saving pictures
            if save:
                plt.savefig(f"frames\gif_{t:05d}.png", dpi=200)

        plt.show()
