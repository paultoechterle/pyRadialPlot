from matplotlib.axes import Axes
from matplotlib.projections import register_projection
from matplotlib.patches import Arc
from matplotlib import collections  as mc
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import MaxNLocator

LAMBDA = 1.55125e-10
G = 0.5


class Radialplot(Axes):

    """
    A RadialPlot or Galbraith Plot
    """

    name = "radialplot"

    LAMBDA = 1.55125e-4
    ZETA = 350
    RHOD = 1.304

    def radialplot(self, Ns, Ni, zeta, rhod, 
                   Dpars=None, marker="o", 
                   transform="Logarithmic"):
       
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
            
        im=self.scatter(self.x, self.y, c=Dpars, marker=marker)
        if Dpars:
            self.figure.colorbar(im, ax=self, orientation="horizontal")
        
        self._add_radial_axis()
        self._add_values_indicators()

    def set_xticks(self, ticks=None):
        if ticks:
            super(Radialplot, self).set_xticks(ticks)
        else:
            loc = LinearLocator(5)
            ticks = loc.tick_values(0., self.max_x)
            super(Radialplot, self).set_xticks(ticks)
        self.spines["bottom"].set_bounds(ticks[0], ticks[-1])

    def set_xlim(self, xlim=None):
        if xlim:
            super(Radialplot, self).set_xlim(xlim[0], 1.25 * xlim[-1])
        else:   
            super(Radialplot, self).set_xlim(0, 1.25 * self.max_x)

    @property
    def max_x(self):
        return np.max(self.x)
        
    @property
    def z(self):
        """ Return transformed z-values"""
        if self.transform == "Linear":
            return  1.0 / LAMBDA * np.log(1.0 + G * self.zeta * LAMBDA * self.rhod * (self.Ns / self.Ni))

        if self.transform == "Logarithmic":
            return np.log(G * self.zeta * LAMBDA * self.rhod * (self.Ns / self.Ni))
           
        if self.transform == "arcsine":
            return np.arcsin(np.sqrt((self.Ns + 3.0/8.0) / (self.Ns + self.Ni + 3.0 / 4.0)))
        
    @property
    def sez(self):
        """Return standard errors"""
        
        if self.transform == "Linear":
            return self.z * np.sqrt( 1.0 / self.Ns + 1.0 / self.Ni)

        if self.transform == "Logarithmic":
            return np.sqrt(1.0 / self.Ns + 1.0 / self.Ni)

        if self.transform == "arcsine":
            return 1.0 / (2.0 * np.sqrt(self.Ns + self.Ni))
        
    @property
    def z0(self):
        """ Return central age"""
        
        if self.transform == "Linear":
            return np.sum(self.z / self.sez**2) / np.sum(1 / self.sez**2)

        if self.transform == "Logarithmic":
            totalNs = np.sum(self.Ns)
            totalNi = np.sum(self.Ni)
            return np.log(G * self.zeta * LAMBDA * self.rhod * (totalNs / totalNi))

        if self.transform == "arcsine":
            return np.arcsin(np.sqrt(np.sum(self.Ns) / np.sum(self.Ns + self.Ni)))
    
    @property
    def x(self):
        return  1.0 / self.sez
    
    @property
    def y(self):
        return (self.z - self.z0) / self.sez
    
    def _z2t(self, z):
        
        if self.transform == "Linear":
            t = z
            return t * 1e-6
        elif self.transform == "Logarithmic":
            NsNi = np.exp(z) / (self.zeta * G * LAMBDA * self.rhod)
        elif self.transform == "arcsine":
            NsNi = np.sin(z)**2 / (1.0 - np.sin(z)**2)
    
        t = 1.0 / LAMBDA * np.log(1.0 + G * self.zeta * LAMBDA * self.rhod * (NsNi))
        return t * 1e-6
    
    def _t2z(self, t):
        
        if self.transform == "Linear":
            return t
        elif self.transform == "Logarithmic":
            return np.log(np.exp(LAMBDA * t) - 1)
        elif self.transform == "arcsine":
            return np.arcsin(
                    1.0 / np.sqrt(
                        1.0 + LAMBDA * self.zeta * G * self.rhod / (np.exp(LAMBDA * t) - 1.0)
                        )
                    )
    
    def _add_radial_axis(self):
        # Get min and max angle
        zr = self._get_radial_ticks_z()

        theta1 = np.arctan(np.min(zr))
        theta1 = np.rad2deg(theta1)
        theta2 = np.arctan(np.max(zr))
        theta2 = np.rad2deg(theta2)

        print(theta1, theta2)
        width = 2.0 * (1.2 * self.max_x)
        height = width

        # The circle is always centered around 0.
        # Width and height are equals (circle)
        arc_element = Arc(
            (0, 0), width, height, theta1=theta1,
            theta2=theta2, linewidth=1, zorder=0, color="k")

        self.add_patch(arc_element)
        
        # Add ticks
        self._add_radial_ticks()
        self._add_radial_ticks_labels()
        
    def _add_radial_ticks(self, nticks=10):

        zr = self._get_radial_ticks_z()

        # Lets build a line collection
        R1 = 1.2 * self.max_x
        R2 = 1.01 * R1
        zr = np.arctan(zr)
        x1 = R1 * np.cos(zr)
        y1 = R1 * np.sin(zr)
        x2 = R2 * np.cos(zr)
        y2 = R2 * np.sin(zr)
        
        starts = list(zip(x1, y1))
        ends = list(zip(x2, y2))
        segments = zip(starts, ends)

        lc = mc.LineCollection(segments, colors='k', linewidths=1)
        self.add_collection(lc)
        
    def _get_radial_ticks_z(self):
        # Let's build the ticks of the Age axis
        za = self.set_radial_ticks_ages()
        zr = self._t2z(np.array(za) * 1e6) - self.z0
        return zr
    
    def set_radial_ticks_ages(self, ticks=None):
        if not ticks:
            ages = self._z2t(self.z)
            start, end = np.int(np.rint(min(ages))), np.int(np.rint(max(ages)))
            loc = MaxNLocator()
            ticks = loc.tick_values(start, end)
        return ticks
        
    def _add_values_indicators(self):
        R1 = (1.2 - 0.02) * self.max_x
        R2 = (1.2 - 0.01) * self.max_x
        ratio = np.arctan(self.y / self.x)
        x1 = R1 * np.cos(ratio)
        y1 = R1 * np.sin(ratio)
        x2 = R2 * np.cos(ratio)
        y2 = R2 * np.sin(ratio)

        starts = list(zip(x1, y1))
        ends = list(zip(x2, y2))
        segments = zip(starts, ends)

        lc = mc.LineCollection(segments, colors='k', linewidths=2)
        self.add_collection(lc) 
        
    def _add_radial_ticks_labels(self):
        # text label
        R3 = 1.2 * self.max_x
        R3 += 0.02 * 1.2 * self.max_x
        za = self.set_radial_ticks_ages()
        labels = self._t2z(np.array(za) * 1e6)
        labels -= self.z0
        x1 = R3 * np.cos(np.arctan(labels))
        y1 = R3 * np.sin(np.arctan(labels))

        for idx, val in enumerate(za):
            self.text(x1[idx], y1[idx], str(val)+ "Ma") 
            
register_projection(Radialplot)


def radialplot(Ns, Ni, zeta, rhod, Dpars=None, marker="o", transform="Logarithmic"):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection="radialplot")
    ax.radialplot(Ns, Ni, zeta, rhod, Dpars, transform=transform)
    return ax
