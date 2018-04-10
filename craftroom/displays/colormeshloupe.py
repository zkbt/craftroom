from .loupe import *

class colormeshloupe(loupe):
    def setup(self, image,
                    ok=None,
                    xaxis=None, yaxis=None,
                    title='', # title to go over the top of the plot
                    figsize=(8,4), # kwargs for figure,
                    dpi=100,
                    hspace=0, wspace=0, # kwargs for gridspec
                    left=0.05, right=0.95,
                    bottom=0.05, top=0.95,
                    width_ratios=[1.0, 0.1],
                    height_ratios=[0.1, 1.0],
                    initialcrosshairs=[0.0, 0.0], # where the crosshairs should start
                    aspect='equal', # kwargs for imshow
                    vmin=None, vmax=None, scale='symlog',
                    labelfontsize=5,
                    datacolor='darkorange',
                    crosshaircolor='darkorange',
                    **kwargs):
        '''
        Initialize the loupe
        (this has a pretty big overhead,
        so use "updateImage" if you can to update
        the data being displayed)
        '''

        self.ok = ok
        self.image = image
        # set the image
        if xaxis is not None:
            self.xaxis = xaxis
        else:
            xsize = self.imagetoplot.shape[1]
            self.xaxis = np.arange(xsize)
        if yaxis is not None:
            self.yaxis = yaxis
        else:
            ysize = self.imagetoplot.shape[0]
            self.yaxis = np.arange(ysize)

        self.dx = np.median(np.diff(self.xaxis))
        self.dy = np.median(np.diff(self.yaxis))

        self.speak('setting up loupe')


        # keep track of where the crosshair is pointing
        self.crosshair = dict(x=initialcrosshairs[0], y=initialcrosshairs[1])

        # set up the figure
        plt.ion()
        self.figure = plt.figure(figsize=figsize, dpi=dpi)
        # create an interactive plot
        self.iplot = iplot(2,2, hspace=hspace, wspace=wspace,
                              left=left, right=right,
                              bottom=bottom, top=top,
                              width_ratios=width_ratios,
                              height_ratios=height_ratios)

        # a dictionary to store the axes objects
        self.ax = {}

        # for displaying the 2D image
        labelkw = dict(fontsize = labelfontsize)
        self.ax['2d'] = self.iplot.subplot(1,0)
        plt.setp(self.ax['2d'].get_xticklabels(), **labelkw)
        plt.setp(self.ax['2d'].get_yticklabels(), **labelkw)

        # for displaying cuts along the dispersion direction
        self.ax['slicex'] = self.iplot.subplot(0,0,sharex=self.ax['2d'])
        self.ax['slicex'].set_title(title, fontsize=8)
        plt.setp(self.ax['slicex'].get_xticklabels(), visible=False)
        plt.setp(self.ax['slicex'].get_yticklabels(), **labelkw)

        # for display cuts along the cross-dispersion direction
        self.ax['slicey'] = self.iplot.subplot(1,1,sharey=self.ax['2d'])
        self.ax['slicey'].xaxis.tick_top()
        self.ax['slicey'].xaxis.set_label_position('top')
        plt.setp(self.ax['slicey'].get_xticklabels(), rotation=270, **labelkw)
        plt.setp(self.ax['slicey'].get_yticklabels(), visible=False)

        # set the limits of the plot
        ok = np.isfinite(self.imagetoplot)
        #self.vmin, self.vmax = np.percentile(self.imagetoplot[ok], [0,100])

        # pick a scale for the plotting
        if scale=='symlog':
            norm = colors.SymLogNorm(
                                  linthresh=craftroom.oned.mad(self.imagetoplot),
                                  linscale=0.1,
                                  vmin=vmin,
                                  vmax=vmax)
        else:
            norm=None


        self.plotted = {}
        # plot the image
        self.plotted['2d'] = self.ax['2d'].pcolormesh(self.xaxis, self.yaxis,
                                                    self.imagetoplot,
                                                  cmap='gray',
                                                  shading='flat',
                                                  zorder=0,
                                                  norm=None)


        # set the x and y limits
        #self.ax['2d'].set_xlim(self.extent[0:2])
        #self.ax['2d'].set_ylim(self.extent[2:4])

        # add crosshair
        crosskw = dict(alpha=0.5, color=crosshaircolor, linewidth=1)
        self.plotted['crossy'] = self.ax['2d'].axvline(self.crosshair['x'],
                                                        **crosskw)
        self.plotted['crossyextend'] = self.ax['slicex'].axvline(
                                                        self.crosshair['x'],
                                                        linestyle='--',
                                                        **crosskw)
        self.plotted['crossx'] = self.ax['2d'].axhline(self.crosshair['y'],
                                                        **crosskw)
        self.plotted['crossxextend'] = self.ax['slicey'].axhline(
                                                        self.crosshair['y'],
                                                        linestyle='--',
                                                        **crosskw)
        # plot slicey
        slicekw = dict(color=datacolor, linewidth=1)

        badalpha = 0.15
        h, v = self.slicey
        ok = self.ok_slicey
        self.plotted['slicey'] = self.ax['slicey'].plot(h[ok], v[ok], **slicekw)[0]
        self.plotted['slicey_bad'] = self.ax['slicey'].plot(h[ok==False], v[ok==False], alpha=badalpha, **slicekw)[0]

        h, v = self.slicex
        ok = self.ok_slicex
        self.plotted['slicex'] = self.ax['slicex'].plot(h[ok], v[ok], **slicekw)[0]
        self.plotted['slicex_bad'] = self.ax['slicex'].plot(h[ok==False], v[ok==False], alpha=badalpha, **slicekw)[0]

        # set the limits of the color scale and the plots
        for a in self.ax.values():
            a.set_autoscaley_on(False)
        self.set_limits(vmin, vmax)

def test():
    fake = createTestImage()
    l = colormeshloupe()
    l.setup(fake)
    l.run()
    return fake, l

def testMovie():
    fake = createTestImage()
    l = colormeshloupe()
    l.setup(fake)
    l.set_limits(0, np.max(fake))
    l.movieSlice(stride=1)
    return fake, l
