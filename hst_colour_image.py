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

def make_image():
    """ Produces HST color image. """
    wdir = "/home/kadu/Dropbox/etacarinae"
    f333w_file = os.path.join(wdir, "j8ma7a0b0_drz.fits")
    f555m_file = os.path.join(wdir, "j8ma7a0e0_drz.fits")
    # Calibrating images
    mags = []
    for imfile in [f333w_file, f555m_file]:
        data = fits.getdata(imfile, 0)
        photflam = fits.getval(imfile, "PHOTFLAM", 1)
        photzpt = fits.getval(imfile, " PHOTZPT", 1)
        mag = -2.5 * np.log10(data * photflam) + photzpt
        mags.append(mag)
    color = mags[1] - mags[0]
    fig = plt.figure(1, figsize=(6, 5.5))
    ax = plt.subplot(111, aspect="equal")
    ax.minorticks_on()
    ax.set_ylim(40, 900)
    ax.set_xlim(50, 800)
    ax.set_xlabel("X (pix)")
    ax.set_ylabel("Y (pix)")
    im = ax.imshow(color, vmin=-.6, vmax=.6, cmap="Spectral_r",
                   origin="bottom")
    cbar = plt.colorbar(im)
    cbar.set_label("F330W - F550M")
    plt.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.98)
    plt.savefig(os.path.join(wdir, "hst_f330w_f550m.png"), dpi=300)

if __name__ == "__main__":
    make_image()