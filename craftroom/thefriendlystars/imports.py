import warnings

import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np


from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u, astropy.coordinates as coord
from astropy.time import Time
from astropy.table import Table


import astroquery.gaia
import astroquery.mast

from tqdm import tqdm


from ..Talker import Talker
