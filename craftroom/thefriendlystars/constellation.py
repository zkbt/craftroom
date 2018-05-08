from .imports import *

# a shortcut getting the coordinates for an object, by its name
get = coord.SkyCoord.from_name


class Constellation(Talker):
    '''
    A Constellation is collection of stars
    that can be accessed through a table of
    astropy coordinates, and/or plotted
    on a Finder chart Panel.
    '''

    def __init__(self,  center,
                        radius=3*u.arcmin,
                        magnitudelimit=20,
                        **kw):
        '''
        Initialize at Constellation object.

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

        # set up the talker for this catalog
        Talker.__init__(self)

        # define the geomery of this Constellation
        self.center = center
        self.radius = radius

        # populate this catalog from its archive source
        self.coneSearch(center=self.center,
                        radius=radius,
                        magnitudelimit=magnitudelimit,
                        **kw)

        # summarize the stars in this constellation
        self.speak('{} contains {} objects'.format(self.name, len(self.objects)))

        # set up a standardized table (a KLUDGE, maybe this should happen earlier?)
        magnames, magdata = [], []
        for k in self.magnitudes.keys():
            magnames.append(k)
            magdata.append(self.magnitudes[k])
        self.table = Table([self.objects]+magdata, names=['coordinates']+magnames)
        self.table['separation'] = self.objects.separation(self.center).arcsec

    def atEpoch(self, epoch=2000):
        '''
        Return SkyCoords of the objects, propagated to a given epoch.

        Parameters
        ----------
        epoch : Time, or float
            Either an astropy time, or a decimal year of the desired epoch.

        Returns
        -------
        coordinates : SkyCoord(s)
            The coordinates, propagated to the given epoch,
            with that epoch stored in the obstime attribute.
        '''

        # calculate the time offset from the epochs of the orignal coordinates
        try:
            # if epoch is already a Time object, pull its decimalyear
            year = epoch.decimalyear
        except AttributeError:
            # if anyting else, assume it is a decimalyear (float or string)
            year = epoch
        newobstime = Time(year, format='decimalyear')
        dt = newobstime - self.objects.obstime

        # calculate the new positions, propagated linearly by dt
        try:
            # if proper motions exist
            newra = self.objects.ra + self.objects.pm_ra_cosdec/np.cos(self.objects.dec)*dt
            newdec = self.objects.dec + self.objects.pm_dec/np.cos(self.objects.dec)*dt
        except TypeError:
            # assume no proper motions, if they're not defined
            newra = self.objects.ra
            newdec = self.objects.dec
            self.speak('no proper motions were used for {}'.format(self.name))

        # return as SkyCoord object
        return coord.SkyCoord(ra=newra, dec=newdec, obstime=newobstime)

    def plot(self, epoch=2000.0, sizescale=10, color=None, **kw):
        '''
        Plot the ra and dec of the coordinates,
        at a given epoch, scaled by their magnitude.

        Parameters
        ----------
        epoch : Time, or float
            Either an astropy time, or a decimal year of the desired epoch.
        sizescale : (optional) float
            The marker size for scatter for a star at the magnitudelimit.
        color : (optional) any valid color
            The color to plot (but there is a default for this catalog.)
        **kw : dict
            Additional keywords will be passed on to plt.scatter.

        Returns
        -------

        plotted : outputs f
        '''

        # pull out the coordinates at this epoch
        coords = self.atEpoch(epoch)

        # calculate the sizes of the stars (logarithmic with brightness?)
        size = np.maximum(sizescale*(1 + self._magnitudelimit - self.magnitude), 0)

        # make a scatter plot of the RA + Dec
        return plt.scatter(coords.ra, coords.dec,
                                    s=size,
                                    color=color or self.color,
                                    label=self.name,
                                    **kw)

    def crossMatchTo(self, reference, radius=1*u.arcsec, visualize=False):
        '''
        Cross-match this catalog onto another reference catalog.
        If proper motions are included in the reference, then
        its coordinates will be propagated to the obstime/epoch
        of this current catalog.

        Parameters
        ----------

        reference : Constellation
            A reference Constellation to which we want to
            cross-match the stars in this catalog. Most likely,
            you'll want to use Gaia for this (since it has
            good astrometry and good proper motions).

        radius : float, with astropy units of angle
            How close to objects need to be in the cross-match
            for them to be considered a matched pair?

        Returns
        -------

        i_this : array of indices
            The elements of this catalog that are matched.

        i_ref : array of indices
            The elements of the reference catalog, corresponding to
        '''

        # find the closest match for each of star in this constellation
        i_ref, d2d_ref, d3d_ref = self.objects.match_to_catalog_sky(reference.atEpoch(self.objects.obstime))

        # extract only those within the specified radius
        ok = d2d_ref < radius
        self.speak('found {} matches within {}'.format(np.sum(ok), radius))

        # make a plot, if desired
        if visualize:
            self.speak('p')
            plt.hist(d2d_ref.arcsec, range=(0,15))
            plt.axvline(radius.arcsec)
            plt.xlabel('Separation (arcsec)')
            plt.ylabel('Number of Matched Sources')

        # return the indices (of this, and of the reference) for the matches
        return ok, i_ref[ok]


class Gaia(Constellation):
    '''
    Gaia catalog contains sources from Gaia DR2,
    including proper motions and parallaxes.
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

        # what's the limiting magnitude for this search?
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


class GALEX(Constellation):
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

class TIC(Constellation):
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
