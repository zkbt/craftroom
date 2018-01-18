'''Tools for resampling array from grid of independent variables to another.'''

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

def binsizes(x):
	'''If x is an array of bin centers, calculate what their sizes are.
		(assumes outermost bins are same size as their neighbors)'''

	binsize = np.zeros_like(x)
	binsize[0:-1] = (x[1:] - x[0:-1])
	binsize[-1] = binsize[-2]
	return binsize

def plotboxy(x, y, **kwargs):
	'''
	Plot with boxes, to show the left and right edges of a box. This is useful,
	for example, to plot flux associated with pixels, in case you are trying to
	do a sub-pixel resample or interpolation or shift.

	(kwargs are passed on to plt.plot)
	'''

	# what are the edges of the bins (making a guess for those on the ends)
	xbinsize = binsizes(x)
	xleft = x - xbinsize/2.0
	xright = x + xbinsize/2.0

	# create a array that doubles up the y values, and interleaves the edges
	plot_x = np.vstack((xleft,xright)).reshape((-1,),order='F')
	plot_y = np.vstack((y,y)).reshape((-1,),order='F')

	# plot those constructed arrays
	plt.plot(plot_x, plot_y, **kwargs)

def fluxconservingresample(xin, yin, xout, test=False, visualize=False, demo=False,
							treatnanas=0.0):
	'''
	Starting from some initial x and y, resample onto a different grid
	(either higher or lower resolution), while conserving total flux.

	When including the entire range of xin, sum(yout) == sum(yin) should be true.
	'''

	# set up the bins, to calculate cumulative distribution of y?
	xinbinsize = binsizes(xin)
	xinleft = xin - xinbinsize/2.0
	xinright = xin + xinbinsize/2.0

	# the first element should be the left edge of the first pixel
	# last element will be right edge of last pixel
	xinforcdf = np.hstack([xinleft, xinright[-1]])

	# to the left of the first pixel, assume flux is zero
	yinforcdf = np.hstack([0, yin])

	# correct for non-finite
	bad = np.isnan(yinforcdf)
	yinforcdf[bad] = treatnanas

	# calculate the cumulative distribution function of the flux (at pixel edge locations)
	cdfin = np.cumsum(yinforcdf)

	# create an interpolator for that
	cdfinterpolator = scipy.interpolate.interp1d(xinforcdf, cdfin,
						kind='linear',
						bounds_error=False,
						fill_value=(0.0, np.sum(yin)))

	# calculate bin edges (of size len(xout)+1)
	xoutbinsize = binsizes(xout)
	xoutleft = xout - xoutbinsize/2.0
	xoutright = xout + xoutbinsize/2.0
	xoutcdf = np.hstack([xoutleft, xoutright[-1]])

	# interpolate the CDF onto those bin edges
	cdfout = cdfinterpolator(xoutcdf)

	# take the derivative of the CDF to get the flux per resampled bin
	# (xout is the center of the bin, and yout is the flux in that bin)
	yout = np.diff(cdfout)

	if visualize:
		fi, (ax_cdf, ax_pdf) = plt.subplots(2,1, sharex=True, figsize=(9,6))
		inkw = dict(color='black', alpha=1, linewidth=3, marker='.', markeredgecolor='none')
		outkw = dict(color='darkorange', alpha=1, linewidth=1, marker='.', markersize=8, markeredgecolor='none')


		legkw = dict(fontsize=10, frameon=False, loc='upper left', bbox_to_anchor=(1, 1))

		# plot the PDFs
		plt.sca(ax_pdf)
		plt.ylabel('Flux per (Original) Pixel')
		plt.xlabel('Pixel')
		# plot the original pixels (in df/dpixel to compare with resampled)
		plotboxy(xin, yin/xinbinsize,
					label='Original Pixels', **inkw)


		# what would a bad interpolation look like?
		badinterpolation = scipy.interpolate.interp1d(xin, yin/xinbinsize,
							kind='linear', bounds_error=False, fill_value=0.0)
		plt.plot(xout, badinterpolation(xout),
						color='cornflowerblue', alpha=1,
						linewidth=1, marker='.',
						markersize=8, markeredgecolor='none',
						label='Silly Simple Interpolation')

		# plot the flux-conserving resampled data (again, in df/d"pixel")
		plt.plot(xout, yout/xoutbinsize, label='Flux-Conserving Interpolation', **outkw)
		plt.legend(**legkw)

		# plot the CDFs
		plt.sca(ax_cdf)
		plt.ylabel('Cumulative Flux (from left)')

		# plot the original CDF
		plt.plot(xinforcdf, cdfin,
					label='Original Pixels', **inkw)

		# plot the interpolated CDF
		plt.plot(xoutcdf, cdfout,
			label='Flux-Conserved Resample', **outkw)
		#plt.legend(**legkw)
		if demo:
			a = raw_input("Pausing a moment to check on interpolation; press return to continue.")

		plt.tight_layout(rect=[0.0, 0.0, 0.67, 1])
		print('{:>6} = {:.5f}'.format('Actual', np.sum(yin)))
		print('{:>6} = {:.5f}'.format('Silly', np.sum(badinterpolation(xout)*xoutbinsize)))
		print('{:>6} = {:.5f}'.format('CDF', np.sum(yout)))


	# return the resampled y-values
	return yout

def testFCR(supersample=True):
	'''this function tests out the resampling code

			supersample=True
				means that there will be multiple new pixels per original pixel

			supersample=False
				means that there will be fewer new pixels than original pixels
	'''

	xinitial = np.arange(39,47)
	yinitial = np.random.uniform(0.0, 0.1, len(xinitial))
	if supersample:
		xresample = np.linspace(np.min(xinitial) - 1.5, np.max(xinitial) + 1.5)
	else:
		xresample = np.linspace(-1,8,5)
	yresample = fluxconservingresample(xinitial, yinitial, xresample, visualize=True)
	plt.savefig('interpolationdemo.pdf')
