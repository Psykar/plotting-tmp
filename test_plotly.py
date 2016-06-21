import plotly.offline as py

import random

import numpy as np
import plotly.graph_objs as go

import scipy.constants
import scipy.integrate
from math import sqrt
from astropy.cosmology import FlatLambdaCDM
from astropy import units
from math import floor
from astropy.coordinates import SkyCoord


class LightCone():
    def __init__(self, box_size=50, ra_min=0, ra_max=np.pi / 2, dec_min=0, dec_max=np.pi / 2, min_z=0, max_z=1, seed=None, hubble=71.0, omega_m=0.27, omega_d=0.73, unique=False, subcone_id=0):
        box_size = box_size
        self.box_size = box_size
        self.ra_min = ra_min * units.rad
        self.ra_max = ra_max * units.rad
        self.dec_min = dec_min * units.rad
        self.dec_max = dec_max * units.rad
        self.min_z = min_z
        self.max_z = max_z
        if seed is None:
            seed = random.randint(0, 1000000)
        self.seed = seed
        self.unique = unique

        self.cosmology = FlatLambdaCDM(hubble, omega_m, omega_d)
        self.subcone_id = subcone_id

    @property
    def min_dist(self):
        retval = self.cosmology.comoving_distance(self.min_z)
        # return redshift_to_comoving_distance(self.min_z)
        return retval

    @property
    def max_dist(self):
        return self.cosmology.comoving_distance(self.max_z)

    @property
    def dec(self):
        return self.dec_max - self.dec_min

    @property
    def ra(self):return self.ra_max - self.ra_min

    @property
    def dec_offset(self):
        return - self.dec_min

    @property
    def ra_offset(self):
        if not self.unique:
            return 0
        box_size = units.Mpc * self.box_size
        comparison = self.max_dist - self.min_dist * np.cos(self.ra)
        if (comparison <= box_size):
            return -self.ra_min

        def angle_root(x):
            phi = self.ra + x
            multiplier = (np.cos(x) - np.sin(x) / np.tan(phi))
            result = box_size - self.max_dist * multiplier

            single = 1 * units.Mpc
            return float(result / single)

        try:
            result = scipy.optimize.ridder(angle_root, 0.5 * np.pi, 0.0)
        except ValueError:
            return None

        return result

    @property
    def origin(self):
        if not self.unique:
            return 0, 0, 0

        theta = self.ra_offset
        phi = theta + self.ra

        # calculate the cone RA height and declination height
        h = self.max_dist * np.sin(phi) - self.min_dist * np.sin(theta)
        h_dec = self.max_dist * np.sin(self.dec)

        # How many will fit in the domain?
        ny = floor(self.box_size / h)

        return (
            self.min_dist * np.cos(phi),
            (self.subcone_id % ny) * h - self.min_dist * np.sin(theta),
            (self.subcone_id / ny) * h_dec
        )

    def tiles(self):
        #  Locate the first box to use. Start by finding a close corner
        #  of the lightcone.
        min_ra = self.ra_min + self.ra_offset
        min_dec = self.dec_min + self.dec_offset

        coord_tmp = SkyCoord(ra=min_ra, dec=min_dec, distance=self.min_dist)
        coord_tmp = coord_tmp.cartesian.x, coord_tmp.cartesian.y, coord_tmp.cartesian.z
        coords_tmp = (sum(x) for x in zip(self.origin, coord_tmp))

        # Make the box corner line up.
        # First is the box number?
        # Coords should always be within a single box
        first = []
        coords = []
        for coord in coords_tmp:
            if coord >= 0:
                first.append(coord / self.box_size)
                coords.append(coord)
            else:

                first.append((coord / self.box_size) - 1)
                coords.append(first[-1] * self.box_size)

        # The goal is to get coords that are -self.box_size off overlapping a box
        ok = self.overlap(coords, (x + self.box_size for x in coords))

    # def overlap(self, min, max):


    def plot(self, ra, dec):
        # Convert ra / dec to radians
        ra = ra * np.pi / 180
        dec = dec * np.pi / 180

        origx, origy, origz = self.origin

        accuracy = 360
        h = self.box_size
        # Calculate cone size from min / max dec, ignore ra?
        cone_angle = (self.dec_max - self.dec_min)
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

        # surface = go.Surface(x=x, y=y, z=z)
        # data = [surface]
        #
        # layout = go.Layout(
        #     title='Parametric Plot',
        #     scene=dict(
        #         xaxis=dict(
        #             gridcolor='rgb(255, 255, 255)',
        #             zerolinecolor='rgb(255, 255, 255)',
        #             showbackground=True,
        #             backgroundcolor='rgb(230, 230,230)',
        #             # range=[0, self.box_size],
        #         ),
        #         yaxis=dict(
        #             gridcolor='rgb(255, 255, 255)',
        #             zerolinecolor='rgb(255, 255, 255)',
        #             showbackground=True,
        #             backgroundcolor='rgb(230, 230,230)',
        #             # range=[0, self.box_size],
        #         ),
        #         zaxis=dict(
        #             gridcolor='rgb(255, 255, 255)',
        #             zerolinecolor='rgb(255, 255, 255)',
        #             showbackground=True,
        #             backgroundcolor='rgb(230, 230,230)',
        #             # range=[0, self.box_size],
        #         ),
        #         aspectmode='cube',
        #     )
        # )
        #
        # fig = go.Figure(data=data, layout=layout)

        boxx= [0, 0, 20, 20, 0, 0, 20, 20],
        boxy= [0, 20, 20, 0, 0, 20, 20, 0],
        boxz= [0, 0, 0, 0, 20, 20, 20, 20],
        boxi= [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2],
        boxj= [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
        boxk= [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6],

        py.plot([
            dict(x=x.value, y=y.value, z=z.value, type='surface'),
            dict(x=boxx, y=boxy, z=boxz, i=boxi, j=boxj, k=boxk, type='mesh3d'),
        ], filename="Light_Cone.html")


if __name__ == '__main__':
    cone = LightCone()
    # tiles = cone.tiles()
    cone.plot(45, 45)
