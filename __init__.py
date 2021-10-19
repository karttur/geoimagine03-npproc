"""
Created 29 Mar 2021

modis
==========================================

Package belonging to Karttur´s GeoImagine Framework.

Author
------
Thomas Gumbricht (thomas.gumbricht@karttur.com)

"""
from .version import __version__, VERSION, metadataD

from .npproc import ProcessNumpy

from geoimagine.modis.modispolar import ProcessModisEase2N