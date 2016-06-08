import plotly.offline as py

from plotly.tools import FigureFactory as FF

import numpy as np
from scipy.spatial import Delaunay

class LightCone():
    def __init__(self, box_size=50, ra_min=0, ra_max=10, dec_min=0, dec_max=10):
        self.box_size = box_size
        self.ra_mix = ra_min
        self.ra_max = ra_max
        self.dec_min = dec_min
        self.dec_max = dec_max

    def plot(self):


        u = np.linspace(0, np.pi * 2, 100)
        v = np.linspace(0, np.pi * 2, 100)
        u, v = np.meshgrid(u, v)
        u = u.flatten()
        v = v.flatten()

        x = u
        y = u * np.cos(v)
        z = u * np.sin(v)

        points2D = np.vstack([u, v]).T
        tri = Delaunay(points2D)
        simplices = tri.simplices

        # define a function for the color assignment
        def dist_from_x_axis(x, y, z):
            return abs(x)

        fig1 = FF.create_trisurf(x=x, y=y, z=z,
                                 colormap=['rgb(255, 10, 120)', 'rgb(255, 100, 255)', ],
                                 color_func=dist_from_x_axis,
                                 simplices=simplices, title="Light Cone",
                                 showbackground=False, gridcolor='rgb(255, 20, 160)',
                                 plot_edges=False, aspectratio=dict(x=1, y=1, z=0.75))
        py.plot(fig1, filename="Light Cone.html")


if __name__ == '__main__':
    cone = LightCone()
    cone.plot()


