from matplotlib.axes import Axes
from matplotlib.projections import register_projection
from matplotlib.patches import Arc
from matplotlib import collections  as mc
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import MaxNLocator
from .radialplot import ZAxis, Radialplot

LAMBDA = 1.55125e-10
G = 0.5

class ZAxisFT(ZAxis):
    
    def _add_radial_axis(self):
        # Get min and max angle

        theta1 = self.ax._t2axis_angle(self.zlim[0] * 1e6)
        theta2 = self.ax._t2axis_angle(self.zlim[1] * 1e6)

        # The circle is always centered around 0.
        # Width and height are equals (circle)
        # Here the easiest is probably to use axis coordinates. The Arc
        # is always centered at (0.,0.) and 

        height = width = 2.0 * self.radius
        arc_element = Arc(
            (0, 0.5), width, height, angle=0., theta1=theta1,
            theta2=theta2, linewidth=1, zorder=0, color="k",
            transform=self.ax.transAxes)

        self.ax.add_patch(arc_element)
        
        # Add ticks
        self.ticks()
        self.labels()
        self.set_zlabel("Age Estimates (Myr)")
        self.add_values_indicators()
    
    def _get_radial_ticks_z(self):
        # Let's build the ticks of the Age axis
        za = self.ticks_locator()
        zr = self.ax._t2z(np.array(za) * 1e6) - self.ax.z0
        return za
    
    def labels(self):
        # text label
        ticks = self.ticks_locator()
        angles = np.array([self.ax._t2axis_angle(val * 1e6) for val in ticks])
        x = 1.02 * self.radius * np.cos(np.deg2rad(angles))
        y = 1.02 * self.radius * np.sin(np.deg2rad(angles)) + 0.5

        for idx, val in enumerate(ticks):
            self.ax.text(x[idx], y[idx], str(val), transform=self.ax.transAxes)

    def ticks(self):

        ticks = self.ticks_locator()
        angles = np.array([self.ax._t2axis_angle(val * 1e6) for val in ticks])
        starts = np.ndarray((len(angles), 2))
        ends = np.ndarray((len(angles), 2))
        starts[:,0] = self.radius * np.cos(np.deg2rad(angles))
        starts[:,1] = self.radius * np.sin(np.deg2rad(angles)) + 0.5
        ends[:,0] = 1.01 * self.radius * np.cos(np.deg2rad(angles))
        ends[:,1] = 1.01 * self.radius * np.sin(np.deg2rad(angles)) + 0.5

        segments = np.stack((starts, ends), axis=1)
        lc = mc.LineCollection(segments, colors='k', linewidths=1, transform=self.ax.transAxes)
        self.ax.add_collection(lc)
    
    def ticks_locator(self, ticks=None):
        if not ticks:
            ages = self.ax._z2t(self.ax.z)
            start, end = np.int(np.rint(min(ages))), np.int(np.rint(max(ages)))
            loc = MaxNLocator()
            ticks = loc.tick_values(start, end)
        return ticks

class FTRadialplot(Radialplot):

    """A RadiaPlot for fission track counts
    
    Returns:
        FTRadialPlot: Radialplot
    """

    name = "fission_track_radialplot"

    LAMBDA = 1.55125e-4
    ZETA = 350
    RHOD = 1.304

    def radialplot(self, Ns, Ni, zeta, rhod, 
                   Dpars=None, transform="logarithmic",
                   **kwargs):
       
        self.Ns = np.array(Ns)
        self.Ni = np.array(Ni)
        Ns = self.Ns[(self.Ns > 0) & (self.Ni > 0)]
        Ni = self.Ni[(self.Ns > 0) & (self.Ni > 0)]
        self.Ns = Ns
        self.Ni = Ni
        self.zeta = zeta
        self.rhod = rhod
        self.Dpars = Dpars
        self.transform = transform
        
        # Prepare the plot Area
        # Left spine
        self.set_ylim(-8, 8)
        self.set_yticks([-2, -1, 0, 1, 2])
        self.spines["left"].set_bounds(-2, 2)
        self.yaxis.set_ticks_position('left')
       
        self.set_xticks()
        self.set_xlim()
        
        self.spines["top"].set_visible(False)
        self.spines["right"].set_visible(False)
            
        im=self.scatter(self.x, self.y, c=Dpars, **kwargs)
        if Dpars:
            self.figure.colorbar(im, ax=self, orientation="horizontal")
        self._add_sigma_lines()
        self._add_central_line()

        self.zaxis = ZAxisFT(self)
        self.zaxis._add_radial_axis()

        # Apply some default labels:
        self.set_xlabel("precision x")
        self.set_ylabel("standardised estimate y")
    
    @property
    def z(self):
        """ Return transformed z-values"""
        if self.transform == "linear":
            return  1.0 / LAMBDA * np.log(1.0 + G * self.zeta * LAMBDA * self.rhod * (self.Ns / self.Ni))

        if self.transform == "logarithmic":
            return np.log(G * self.zeta * LAMBDA * self.rhod * (self.Ns / self.Ni))
           
        if self.transform == "arcsine":
            return np.arcsin(np.sqrt((self.Ns + 3.0/8.0) / (self.Ns + self.Ni + 3.0 / 4.0)))
        
    @property
    def sez(self):
        """Return standard errors"""
        
        if self.transform == "linear":
            return self.z * np.sqrt( 1.0 / self.Ns + 1.0 / self.Ni)

        if self.transform == "logarithmic":
            return np.sqrt(1.0 / self.Ns + 1.0 / self.Ni)

        if self.transform == "arcsine":
            return 1.0 / (2.0 * np.sqrt(self.Ns + self.Ni))
        
    @property
    def z0(self):
        """ Return central age"""
        
        if self.transform == "linear":
            return np.sum(self.z / self.sez**2) / np.sum(1 / self.sez**2)

        if self.transform == "logarithmic":
            totalNs = np.sum(self.Ns)
            totalNi = np.sum(self.Ni)
            return np.log(G * self.zeta * LAMBDA * self.rhod * (totalNs / totalNi))

        if self.transform == "arcsine":
            return np.arcsin(np.sqrt(np.sum(self.Ns) / np.sum(self.Ns + self.Ni)))
    
    def _z2t(self, z):
        
        if self.transform == "linear":
            t = z
            return t * 1e-6
        elif self.transform == "logarithmic":
            NsNi = np.exp(z) / (self.zeta * G * LAMBDA * self.rhod)
        elif self.transform == "arcsine":
            NsNi = np.sin(z)**2 / (1.0 - np.sin(z)**2)
    
        t = 1.0 / LAMBDA * np.log(1.0 + G * self.zeta * LAMBDA * self.rhod * (NsNi))
        return t * 1e-6
    
    def _t2z(self, t):
        
        if self.transform == "linear":
            return t
        elif self.transform == "logarithmic":
            return np.log(np.exp(LAMBDA * t) - 1)
        elif self.transform == "arcsine":
            return np.arcsin(
                    1.0 / np.sqrt(
                        1.0 + LAMBDA * self.zeta * G * self.rhod / (np.exp(LAMBDA * t) - 1.0)
                        )
                    )


register_projection(FTRadialplot)

def radialplot(Ns=None, Ni=None, zeta=None, rhod=None, file=None,
               Dpars=None, transform="logarithmic", **kwargs):
    """Plot Fission Track counts using a RadialPlot (Galbraith Plot)
    
    Args:
        Ns (list or numpy array, optional): 
            Spontaneous counts. 
            Defaults to None.
        Ni (list or numpy array, optional): 
            Induced counts. 
            Defaults to None.
        zeta (float, optional): 
            Zeta calibration parameter.
            Defaults to None.
        rhod (float, optional): 
            Rhod calibration parameter.
            Defaults to None.
        file (string, optional): 
            Data File, for now pyRadialPlot only accepts
            file format similar to RadialPlotter.
            Defaults to None.
        Dpars (list or numpy array or float, optional): 
            Dpars values associated with the grain counts.
            Defaults to None.
        transform (str, optional): 
            Transformation used.
            Options are "linear", "logarithmic", "arcsine".
            Defaults to "logarithmic".
        kwargs: Matplotlib additional parameters.
    
    Returns:
        matplotlib.Axes: 
            A Matplotlib Axes object.
    """
    if not "color" in kwargs.keys():
        kwargs["color"] = "black"

    fig = plt.figure(figsize=(6,6))
    if file:
        from .utilities import read_radialplotter_file
        data = read_radialplotter_file(file)
        Ns = data["Ns"]
        Ni = data["Ni"]
        zeta = data["zeta"]
        rhod = data["rhod"]

    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection="fission_track_radialplot")
    ax.radialplot(Ns, Ni, zeta, rhod, Dpars, transform=transform, **kwargs)
    return ax
