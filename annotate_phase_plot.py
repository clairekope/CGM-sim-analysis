#!/usr/bin/env python
# coding: utf-8

import numpy as np
import yt
import yt.units as u

def annotate_isobars_adiabats(phaseplot, z_field, log_step=2, label_limit_frac=0.95):
    """
    Annotate isobars and adiabats to a yt density-temperature PhasePlot.
    These contours are designed to be plotted every few powers of ten.
    Only a handful of contours are labeled; isobars are labeled along the
    right y-axis and adiabats along the bottom x-axis.
    Nothing is returned; the plot is modified "in place."
    
    phaseplot : the density-temperature PhasePlot
    z_field : name of the z-field (binned field) to be annotated 
              in the PhasePlot
    log_step : control the step between contours; e.g., log_step=2 plots
               every other power of ten
    label_limit_frac : if a contour would be labeled, but the label's anchor
                       falls beyond this fraction of the respective axis,
                       the label is not included. This prevents text 
                       overflowing the axis bounds. Adjust if too many/too
                       few labels are being shown.
    """

    pt = phaseplot.plots[z_field]
    ax = pt.axes

    d_lo, d_hi = ax.get_xlim()
    T_lo, T_hi = ax.get_ylim()

    d_lo *= u.g/u.cm**3
    d_hi *= u.g/u.cm**3
    T_lo *= u.K
    T_hi *= u.K

    logP_lo, logP_hi = np.trunc(
        [np.log10( (d_lo/u.mh * u.kb*T_lo).to('dyn*cm**-2') ),
         np.log10( (d_hi/u.mh * u.kb*T_hi).to('dyn*cm**-2') )]
    )

    logK_lo, logK_hi = np.trunc(
        [np.log10( (u.kb*T_lo/(d_hi/u.mh)**(2/3)).to('eV*cm**2', equivalence='thermal') ),
         np.log10( (u.kb*T_hi/(d_lo/u.mh)**(2/3)).to('eV*cm**2', equivalence='thermal') )]
    )

    d = np.linspace(d_lo, d_hi)
    logd = np.log10(d.value)
    d_label_limit = 10**(logd[0] + label_limit_frac*(logd[-1] - logd[0]))

    logP = np.linspace(logP_lo, logP_hi, int((logP_hi-logP_lo)/log_step)+1)
    logK = np.linspace(logK_lo, logK_hi, int((logK_hi-logK_lo)/log_step)+1)

    for P in 10**logP * u.dyn*u.cm**-2:
        T = P/(u.kb*d/u.mh)
        logT = np.log10(T)
        T_label_limit = 10**(logT[0]+label_limit_frac*(logT[-1] - logT[0]))

        ax.plot(d, T, color='gray', ls=':')

        if T[-1] < T_lo or T[-1] > T_label_limit:
            continue
            
        point1 = ax.transData.transform_point(np.array([d[-2], T[-2]]))
        point2 = ax.transData.transform_point(np.array([d[-1], T[-1]]))
        angle = np.rad2deg(np.arctan2(point2[1]-point1[1], point2[0]-point1[0]))

        ax.text(d[-1], T[-1], f'{P.v:.0e}'+r' dyn cm$^{-2}$   ',
                horizontalalignment='right', verticalalignment='top', fontsize=12,
                rotation=angle, rotation_mode='anchor')

    for K in 10**logK * u.eV*u.cm**2:
        T = (K/u.kb*(d/u.mh)**(2/3)).to('K', equivalence='thermal')
        logT = np.log10(T.value)

        ax.plot(d, T, color='gray', ls=':')

        m = (logT[-1]-logT[-2])/(logd[-1]-logd[-2])
        b = logT[-1] - m*logd[-1]
        x_trcpt = 10**((np.log10(T_lo.v)-b) / m)

        # doesn't intersect x-axis within plot window
        if x_trcpt < d_lo.value or x_trcpt > d_label_limit:
            continue

        point1 = ax.transData.transform_point(np.array([d[-2], T[-2]]))
        point2 = ax.transData.transform_point(np.array([d[-1], T[-1]]))
        angle = np.rad2deg(np.arctan2(point2[1]-point1[1], point2[0]-point1[0]))

        ax.text(x_trcpt, T_lo, f'   {K.v:.0e}'+r' eV cm$^2$',
                horizontalalignment='left', verticalalignment='bottom', fontsize=12,
                rotation=angle, rotation_mode='anchor')
        
