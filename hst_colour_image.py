# -*- coding: utf-8 -*-
"""

Created on 30/10/2017

@Author: Carlos Eduardo Barbosa

Produces color image for Eta Carinae using HST images.

"""
from __future__ import division, print_function

import os

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import ndimage

def make_image():
    """ Produces HST color image. """
    wdir = "/home/kadu/Dropbox/etacarinae"
    f333w_file = os.path.join(wdir, "j8ma7a0b0_drz.fits")
    f555m_file = os.path.join(wdir, "j8ma7a0e0_drz.fits")
    coords = SkyCoord("10h45m03.692s -59d41m04.4",
                      unit=[u.hourangle, u.degree])
    # Calibrating images
    mags = []
    for i, imfile in enumerate([f333w_file, f555m_file]):
        data = fits.getdata(imfile, 0)
        photflam = fits.getval(imfile, "PHOTFLAM", 1)
        photzpt = fits.getval(imfile, " PHOTZPT", 1)
        mag = -2.5 * np.log10(data * photflam) + photzpt
        mags.append(mag)
        header = fits.getheader(imfile, 1)
        wcs = WCS(header)
    x = np.arange(header["NAXIS1"]) + 1
    y = np.arange(header["NAXIS2"]) + 1
    xx, yy = np.meshgrid(x, y)
    xy = np.column_stack((xx.flatten(), yy.flatten()))
    radec = wcs.all_pix2world(xy, 1)
    ra = radec[:,0].reshape(xx.shape) * u.degree
    dec = radec[:,1].reshape(yy.shape) * u.degree
    color = mags[0] - mags[1]
    ra = (ra.value - coords.ra.value) * 3600 / 2
    dec = (dec.value - coords.dec.value) * 3600
    # Starting plot
    plt.style.context("seaborn-paper")
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "serif"
    plt.rcParams['font.serif'] = 'Computer Modern'
    fig = plt.figure(1, figsize=(5,4.1))
    ax = plt.subplot(111)
    ax.minorticks_on()
    ax.tick_params(right=True, top=True, axis="both",
                   direction='in', which="both")
    color[np.isnan(color)] = -9999
    color = ndimage.median_filter(color, 4)
    color[color < -10] = np.nan
    im = ax.pcolormesh(ra, dec,color,
                       vmin=0, vmax=1, cmap="Spectral_r")
    ax.set_aspect("equal")
    radius = 11
    ax.set_ylim(-radius, radius)
    ax.set_xlim(radius, -radius)
    ax.set_xlabel("$\Delta$RA (arcsec)")
    ax.set_ylabel("$\Delta$DEC (arcsec)")
    # ax.invert_yaxis()
    cbar = plt.colorbar(im, fraction=0.048, pad=0.01)
    cbar.set_label("F330W - F550M")
    cbar.ax.tick_params(right=True, top=True, axis="both",
                   direction='in', which="both")
    radius = [0.3, 3, 9.5, 9.95]
    c = ["w", "w", "m", "m"]
    s = ["r=0.15''", "r=3''", "r=9.5''", "r=9.95''"]
    pos = [(-1, 0.5), (2.2,2.2), (1.5,8), (6.2,8)]
    for i,r in enumerate(radius):
        circle = plt.Circle((0,0), radius=r, ec=c[i], fc="none",
                            linestyle="--")
        ax.add_patch(circle)
        # ax.text(-pos[i][0], pos[i][1], s[i], color=c[i], fontsize=8,
        #         fontweight='bold')
    plt.subplots_adjust(left=0.06, right=0.9, bottom=0.1, top=0.98)
    plt.savefig(os.path.join(wdir, "hst_f330w_f550m.png"), dpi=300)
    plt.show()

if __name__ == "__main__":
    make_image()