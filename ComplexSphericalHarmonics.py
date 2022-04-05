#!/usr/bin/env python 3.6.8 Anaconda 64-bit
# -*- coding: utf-8 -*-
#Created on Tue Apr  5 11:44:49 2022
#@author: frankgrijalva
#But I only modified work by @christian from SciPy 15 Jun 2020
#https://scipython.com/blog/visualizing-the-real-forms-of-the-spherical-harmonics/
#
#This is a script to Visualize complex spherical harmonic(s) of your choosing

#Trey's Top-Shelf Imports:
import sys
##########################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# The following import configures Matplotlib for 3D plotting.
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
# I wants this to work, but it refuses to work with LaTeX ∴usetex=False
plt.rc('text', usetex=False)

# Grids of polar and azimuthal angles
theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)
# Create a 2-D meshgrid of (theta, phi) angles.
theta, phi = np.meshgrid(theta, phi)
# Calculate the Cartesian coordinates of each point in the mesh.
xyz = np.array([np.sin(theta) * np.sin(phi),
                np.sin(theta) * np.cos(phi),
                np.cos(theta)])


#Func-y town
def plot_Ysq(ax, el, m):
    """This Function will plot the spherical harmonic of degree el and order m on Axes ax."""

    # NB In SciPy's sph_harm function the azimuthal coordinate, theta,
    # comes before the polar coordinate, phi.
    Ysq = (sph_harm(abs(m), el, phi, theta))**2

#Dr Brookes doesn't want us plotting real QM yet, he wants complex, actual spherical harmonics
    # Linear combination of Y_l,m and Y_l,-m to create the real form.
    # if m < 0:
    #     Y = np.sqrt(2) * (-1)**m * Y.imag
    # elif m > 0:
    #     Y = np.sqrt(2) * (-1)**m * Y.real
    Yx, Yy, Yz = np.abs(Ysq) * xyz

    # Color the plotted surface according to the phase:
    #(While red is expanding blue contracts)
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet'))
    #cmap.set_clim(-0.5, 0.5)

    ax.plot_surface(Yx, Yy, Yz,
                    facecolors=cmap.to_rgba(Ysq.real),
                    rstride=1, cstride=1)

    # Draw a set of x, y, z axes for reference.
    ax_lim = 0.5
    ax.plot([-ax_lim, ax_lim], [0,0], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [-ax_lim, ax_lim], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [0,0], [-ax_lim, ax_lim], c='0.5', lw=1, zorder=10)
    # Set the Axes limits and title, turn off the Axes frame.
    ax.set_title(r'$Y_{{{},{}}}$'.format(el, m))
    ax_lim = 0.25
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.axis('off')

#cleaning any plots leftover from before
plt.close("all")


#catching error input for a specific l&m
try:
     #asking for user input for l and m
     l = int(input("Please enter an int for l, /nthe 'angular momentum quantum number': "))
     m = int(input("Please enter an int for m, /nthe 'magnetic quantum number': "))

     #Plotting:

     ##The 'coolest one', Y_3,0 will be plotted
     ##because everyone loves d and f orbitals
     fig1 = plt.figure(figsize=plt.figaspect(1.))
     ax = fig1.add_subplot(projection='3d')
     plot_Ysq(ax, 3, 0)

     #Now we'll need that user input to plot whichever they wanted to see
     fig2 = plt.figure(figsize=plt.figaspect(1.))
     ax = fig2.add_subplot(projection='3d')
     plot_Ysq(ax, l, m)

     #To plot a family of these functions, uncomment below:
     # el_max = 3
     # figsize_px, DPI = 900, 100
     # figsize_in = figsize_px / DPI
     # fig = plt.figure(figsize=(figsize_in, figsize_in), dpi=DPI)
     # spec = gridspec.GridSpec(ncols=2*el_max+1, nrows=el_max+1, figure=fig)
     # for el in range(el_max+1):
     #     for m_el in range(-el, el+1):
     #         print(el, m_el)
     #         ax = fig.add_subplot(spec[el, m_el+el_max], projection='3d')
     #         plot_Ysq(ax, el, m_el)
     # plt.tight_layout()
     # plt.show()

except: sys.exit("Please run again and enter an appropriate int for l&m")
