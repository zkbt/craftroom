
import matplotlib.pyplot as plt
import numpy as np


from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u, astropy.coordinates as coord
from astropy.time import Time
from astropy.table import Table

from craftroom.star import Star
from craftroom.Talker import Talker
from craftroom.oned import mad


import astroquery.gaia
import astroquery.mast
