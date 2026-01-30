"""
redshift_slider - Interactive Redshift Fitting Tool
====================================================

A simple tool to visually check emission line positions at different redshifts.
Use the slider to move where strong lines should appear on your spectrum.

Basic usage:
    from redshift_slider import RedshiftSlider
    
    slider = RedshiftSlider(wavelength, flux, z_prior=2.5, msaid='12345')
    z_fitted = slider.show()

For more options, see RedshiftSlider docstring.
"""

from .core import RedshiftSlider, STRONG_LINES, ZOOM_LINES, LINE_COLORS

__version__ = "0.1.0"
__author__ = "Raunaq Rai"
__all__ = ["RedshiftSlider", "STRONG_LINES", "ZOOM_LINES", "LINE_COLORS"]
