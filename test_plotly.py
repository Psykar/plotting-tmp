import plotly.offline as py

import random

import numpy as np
import plotly.graph_objs as go


class LightCone():
    def __init__(self, box_size=50, ra_min=0, ra_max=10, dec_min=0, dec_max=10, seed=None):
        self.box_size = box_size
        self.ra_mix = ra_min
        self.ra_max = ra_max
        self.dec_min = dec_min
        self.dec_max = dec_max
        if seed is None:
            seed = random.randint(0, 1000000)
        self.seed = seed

    def calc(self):
        origin = self.calc_origin()


    def calc_origin(self):
        self.calc_subcone_angle()
        self.calc_subcone_origin()

    def calc_subcone_angle(self):
        pass

    def calc_subcone_origin(self):
        pass

    def plot(self, origin, ra, dec):
        # Convert ra / dec to radians
        ra = ra * np.pi / 180
        dec = dec * np.pi / 180

        origx, origy, origz = origin

        accuracy = 360
        h = self.box_size
        # Calculate cone size from min / max dec, ignore ra?
        cone_angle = (self.dec_max - self.dec_min) * np.pi / 180
        r = h * np.sin(cone_angle)

        phi = np.linspace(0, np.pi * 2, accuracy)
        u = np.linspace(0, h, accuracy)
        phigrid, ugrid,  = np.meshgrid(phi, u)

        x = h - ugrid
        y = ((h - ugrid) / h) * r * np.cos(phigrid)
        z = ((h - ugrid) / h) * r * np.sin(phigrid)

        # Dec has to come first otherwise it's rotation wouldn't
        # be around an axis and maths gets complicated
        # Rotate by dec, around y axis
        xnew = x * np.cos(dec) - z * np.sin(dec)
        znew = z * np.cos(dec) + x * np.sin(dec)
        x = xnew
        z = znew

        # Rotate by RA, around z axis
        xnew = x * np.cos(ra) - y * np.sin(ra)
        ynew = y * np.cos(ra) + x * np.sin(ra)
        x = xnew
        y = ynew

        x = x + origx
        y = y + origy
        z = z + origz

        surface = go.Surface(x=x, y=y, z=z)
        data = [surface]

        layout = go.Layout(
            title='Parametric Plot',
            scene=dict(
                xaxis=dict(
                    gridcolor='rgb(255, 255, 255)',
                    zerolinecolor='rgb(255, 255, 255)',
                    showbackground=True,
                    backgroundcolor='rgb(230, 230,230)',
                    range=[0, self.box_size],
                ),
                yaxis=dict(
                    gridcolor='rgb(255, 255, 255)',
                    zerolinecolor='rgb(255, 255, 255)',
                    showbackground=True,
                    backgroundcolor='rgb(230, 230,230)',
                    range=[0, self.box_size],
                ),
                zaxis=dict(
                    gridcolor='rgb(255, 255, 255)',
                    zerolinecolor='rgb(255, 255, 255)',
                    showbackground=True,
                    backgroundcolor='rgb(230, 230,230)',
                    range=[0, self.box_size],
                ),
                aspectmode='cube',
            )
        )

        fig = go.Figure(data=data, layout=layout)
        py.plot(fig, filename="Light_Cone.html")


if __name__ == '__main__':
    cone = LightCone()
    cone.plot((0, 0, 0), 45, 45)


