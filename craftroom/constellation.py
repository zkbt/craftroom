
import matplotlib.pyplot as plt
import numpy as np


from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u, astropy.coordinates as coord
from astropy.time import Time

from craftroom.star import Star
from craftroom.Talker import Talker
from craftroom.oned import mad


import astroquery.skyview
import astroquery.gaia
import astroquery.mast



# define the images that accessible to skyview
twomass = ['2MASS-J', '2MASS-H', '2MASS-K']
ukidss = ['UKIDSS-Y', 'UKIDSS-J', 'UKIDSS-H', 'UKIDSS-K']
wise = ['WISE 3.4', 'WISE 4.6', 'WISE 12', 'WISE 22']
dss1 = ['DSS1 Blue', 'DSS1 Red']
dss2 = ['DSS2 Blue', 'DSS2 Red']
GALEX = ['GALEX Far UV', 'GALEX Near UV']

class Image:
    '''
    This represents images that lines up with a given patch of the sky.
    '''

    def __init__(self, hdu, name=None):
        '''
        Initialize an image.

        Parameters
        ----------

        hdu : a PrimaryHDU file
            FITS file
        '''

        self.header = hdu.header
        self.data = hdu.data
        self.wcs = WCS(hdu.header)
        self.name = name


class Catalog(Talker):
    '''
    A collection of points that can be overplotted.
    '''

    def __init__(self, center, radius=3*u.arcmin, magnitudelimit=20, **kw):
        '''
        Initialize at Catalog object.

        Parameters
        ----------
        center : SkyCoord object
            The center around which the query will be made.
        radius : float, with units of angle
            The angular radius for the query.
        magnitudelimit : float
            The maximum magnitude to include in the download.
            (This is explicitly thinking UV/optical/IR, would
            need to change to flux to be able to include other
            wavelengths.)
        '''

        # set up the talke for this catalog
        Talker.__init__(self)

        # populate this catalog from its archive source
        self.coneSearch(center, radius=radius, magnitudelimit=magnitudelimit, **kw)


        self.speak('{} contains {} objects'.format(self.name, len(self.objects)))



    def atEpoch(self, epoch=2000):
        '''
        Return SkyCoord object propagated to a given epoch.

        Parameters
        ----------

        epoch : Time, or float
            Either an astropy time, or a decimal year of the desired epoch.
        '''

        # calculate the time offset
        try:
            year = epoch.decimalyear
        except AttributeError:
            year = epoch

        newobstime = Time(year, format='decimalyear')
        dt = newobstime - self.objects.obstime

        # calculate the new positions
        try:
            newra = self.objects.ra + self.objects.pm_ra_cosdec/np.cos(self.objects.dec)*dt
            newdec = self.objects.dec + self.objects.pm_dec/np.cos(self.objects.dec)*dt
        except TypeError:
            newra = self.objects.ra
            newdec = self.objects.dec
            self.speak('no proper motions were used for {}'.format(self.name))


        # return as SkyCoord object
        return coord.SkyCoord(ra=newra, dec=newdec, obstime=newobstime)

    def plot(self, epoch=2000.0, sizescale=10, color=None, **kw):

        # pull out the coordinates at this epoch
        coords = self.atEpoch(epoch)

        # calculate the sizes of the stars (logarithmic with brightness?)
        size = -sizescale*(self.magnitude - self._magnitudelimit)
        self.plotted = plt.scatter(coords.ra, coords.dec,
                                    s=size,
                                    color=color or self.color,
                                    label=self.name,
                                    **kw)





class GALEX(Catalog):
    name = 'GALEX'
    color = 'orchid'

    def coneSearch(self, center, radius=3*u.arcmin, magnitudelimit=25):
        '''
        Run a cone search of the GALEX archive
        '''


        self._magnitudelimit = magnitudelimit

        # run the query
        self.speak('querying GALEX, centered on {} with radius {}'.format(center, radius, magnitudelimit))

        coordinatetosearch = '{0.ra.deg} {0.dec.deg}'.format(center)
        self._table = astroquery.mast.Catalogs.query_region(coordinates=center, radius=radius, catalog='GALEX')



        # the gaia DR2 epoch is 2015.5
        epoch = 2005#???

        # create skycoord objects
        self.objects = coord.SkyCoord(  ra=self._table['ra'].data*u.deg,
                                        dec=self._table['dec'].data*u.deg,
                                        obstime=Time(epoch, format='decimalyear'))

        self.magnitudes = dict(NUV=self._table['nuv_mag'].data, FUV=self._table['fuv_mag'].data)
        self.magnitude = self.magnitudes['NUV']

class TIC(Catalog):
    name = 'TIC'
    color = 'green'

    def coneSearch(self, center, radius=3*u.arcmin, magnitudelimit=25):
        '''
        Run a cone search of the GALEX archive
        '''


        self._magnitudelimit = magnitudelimit

        # run the query
        self.speak('querying TIC, centered on {} with radius {}'.format(center, radius, magnitudelimit))

        coordinatetosearch = '{0.ra.deg} {0.dec.deg}'.format(center)
        self._table = Catalogs.query_region(coordinates=center, radius=radius, catalog='TIC')



        # the gaia DR2 epoch is 2015.5
        epoch = 2000#???

        # create skycoord objects
        self.objects = coord.SkyCoord(  ra=self._table['ra'].data*u.deg,
                                        dec=self._table['dec'].data*u.deg,
                                        obstime=Time(epoch, format='decimalyear'))

        self.magnitudes = dict(T=self._table['Tmag'].data)
        self.magnitude = self.magnitudes['T']




class Gaia(Catalog):
    '''
    Create a catalog of Gaia sources, as astropy SkyCoords.
    '''

    name = 'Gaia'
    color = 'black'

    def query(self, query):
        '''
        Send an ADQL query to the Gaia archive,
        wait for a response,
        and hang on to the results.
        '''

        # send the query to the Gaia archive
        self._gaia_job = astroquery.gaia.Gaia.launch_job(query)

        # return the table of results
        return self._gaia_job.get_results()

    def coneSearch(self, center, radius=3*u.arcmin, magnitudelimit=20):
        '''
        Run a cone search of the Gaia DR2.
        '''


        self._magnitudelimit = magnitudelimit

        # define a query for cone search surrounding this center
        conequery = """SELECT source_id,ra,ra_error,dec,dec_error,pmra, pmra_error, pmdec, pmdec_error, parallax, parallax_error, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, radial_velocity, radial_velocity_error, phot_variable_flag, teff_val, a_g_val FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{},{},{}))=1 and phot_g_mean_mag < {}""".format(center.ra.deg, center.dec.deg, radius.to(u.deg).value, magnitudelimit)


        # run the query
        self.speak('querying Gaia DR2, centered on {} with radius {}, for G<{}'.format(center, radius, magnitudelimit))
        self._table = self.query(conequery)

        # tidy up quantities, setting motions to 0 if poorly defined
        for key in ['pmra', 'pmdec', 'parallax', 'radial_velocity']:
            bad = self._table[key].mask
            self._table[key][bad] = 0.0

        # the gaia DR2 epoch is 2015.5
        epoch = 2015.5

        # create skycoord objects
        self.objects = coord.SkyCoord(  ra=self._table['ra'].data*u.deg,
                                        dec=self._table['dec'].data*u.deg,
                                        pm_ra_cosdec=self._table['pmra'].data*u.mas/u.year,
                                        pm_dec=self._table['pmdec'].data*u.mas/u.year,
                                        radial_velocity=self._table['radial_velocity'].data*u.km/u.s,
                                        distance=1000*u.pc/self._table['parallax'].data,
                                        obstime=Time(epoch, format='decimalyear'))

        self.magnitudes = dict(         G=self._table['phot_g_mean_mag'].data,
                                        B=self._table['phot_bp_mean_mag'].data,
                                        R=self._table['phot_rp_mean_mag'].data)

        self.magnitude = self.magnitudes['G']


class TwoMass(Gaia):
    '''A catalog for 2MASS, using the Gaia-archive hosted search.'''

    name = '2MASS - J'
    color = 'orange'
    zorder = -1
    def coneSearch(self, center, radius=3*u.arcmin, magnitudelimit=20):
        '''
        Run a cone search of the Gaia DR2.
        '''


        self._magnitudelimit = magnitudelimit

        # define a query for cone search surrounding this center
        conequery = """SELECT designation, ra, dec, j_m, h_m, ks_m, j_date FROM gaiadr1.tmass_original_valid WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{},{},{}))=1 and j_m < {}""".format(center.ra.deg, center.dec.deg, radius.to(u.deg).value, magnitudelimit)


        # run the query
        self.speak('querying 2MASS, centered on {} with radius {}, for J<{}'.format(center, radius, magnitudelimit))
        self._table = self.query(conequery)




        # create skycoord objects
        self.objects = coord.SkyCoord(  ra=self._table['ra'].data*u.deg,
                                        dec=self._table['dec'].data*u.deg,
                                        obstime=Time(self._table['j_date'].data, format='jd'))

        self.magnitudes = dict(         J=self._table['j_m'].data,
                                        H=self._table['h_m'].data,
                                        Ks=self._table['ks_m'].data)
        self.magnitude = self.magnitudes['J']


# define som
class Constellation(Talker):
    '''
    A constellation object represents a small patch of
    the sky, centered on a particular target.
    '''

    def __init__(self, name, radius=3*u.arcmin):
        self.star = Star(name)
        self.center = self.star.icrs
        self.radius = radius

    def populateImagesFromSurveys(self, surveys=dss2 + twomass):
        '''
        Load images from archives.
        '''

        # what's the coordinate center?
        coordinatetosearch = '{0.ra.deg} {0.dec.deg}'.format(self.center)

        # query sky view for those images
        paths = astroquery.skyview.SkyView.get_images(
                                    position=coordinatetosearch,
                                    radius=self.radius,
                                    survey=surveys)

        # populate the images for each of these
        self.images = [Image(p[0], s) for p, s in zip(paths, surveys)]

    def populateCatalogsFromSurveys(self, surveys=[Gaia, TwoMass, GALEX, TIC]):

        self.catalogs = {}
        for s in surveys:
            this = s(self.center, self.radius)
            self.catalogs[this.name] = this

    def plotGrid(self):

        N = len(self.images)
        fig = plt.figure(figsize=(20, 21*N), dpi=200)
        self.ax = {}
        share = None
        for i, image in enumerate(self.images):
            ax = fig.add_subplot(1, N, i+1, projection=image.wcs)



            norm = plt.matplotlib.colors.SymLogNorm(
                                  linthresh=mad(image.data),
                                  linscale=0.1,
                                  vmin=None,
                                  vmax=None)

            ax.imshow(image.data, origin='lower', cmap='gray_r', norm=norm, alpha=0.5)
            transform = ax.get_transform('world')
            ax.set_title(image.name)
            #ax.grid(color='white', ls='solid')

class Finder:
    def __init__(self, constellation):
        self.constellation = constellation



class Frame:
    '''
    A single frame of a plot,
    that has up to one image in the background,
    and any number of catalogs plotted.
    '''
    def __init__(self, image, catalogs=None):
        pass
        #???
