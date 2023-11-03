import time
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from tqdm import tqdm
import multiprocessing as mp

class Perlin():
    def __init__(self, seed:int=6942021, precision:int=3) -> None:
        self._seed = seed
        # Store jacobian and perlin values to avoid having to recalculate them
        self.perlin_values = {}
        self.jacobian_values = {}
        self.precision = precision
        self.values = []

    def _random_number_generator(self, value1:int, value2:int, value3:int, seed:int) -> int:
        """Random number generator to generate random gradient

        Args:
            value1: x-coordinate contribution to total value
            value2: y-coordinate contribution to total value
            value3: z-coordinate contribution to total value
            seed: seed for random number generation
            
        Returns:
            (int): random number
        """

        total_value = value1 + value2 +  value3
        m_value = np.int64(total_value % 0xFFFFFFFF)
        BIT_NOISE1 = np.int64(0xB5297A4D)
        BIT_NOISE2 = np.int64(0x68E31DA4)
        BIT_NOISE3 = np.int64(0x1B56C4E9)

        m_value += BIT_NOISE1
        m_value += np.int64(seed)
        m_value ^= np.right_shift(m_value, 8, casting='unsafe', dtype='int64')
        m_value += BIT_NOISE2
        m_value ^= np.left_shift(m_value, 8, casting='unsafe', dtype='int64')
        m_value += BIT_NOISE3
        m_value ^= np.right_shift(m_value, 8, casting='unsafe', dtype='int64')
        return int(m_value)

    def _get_lattice_indices(self, x: float, y:float, z:float) -> int:
        """Return the lattice indices of the unit cube containing the given point
                                    (d)-----------(e)
                                    /|            /| 
                                   / |           / |
                                  /  |          /  |
                                 /  (b)--------/--(g)
                               (c)-----------(f)  /
                                |  /          |  /
                                | /           | /
                                |/            |/           
                               (a)-----------(h)
        a = x0, y0, z0      e = x1, y1, z1
        b = x0, y1, z0      f = x1, y0, z1
        c = x0, y0, z1      g = x1, y1, z0
        d = x0, y1, z1      h = x1, y0, z0              
        x0, y0, z0 indicate smaller values; x1, y1, z1 indicate larger values                        
                                
        Args:
            x: float describing x-coordinate
            y: float describing y-coordinate
            z: float describing z-coordinate
        
        Returns:
            x0, y0, z0: integers containing the lower values of the lattice cube
            x1, y1, z1: integers containing the upper values of the lattice cube
        """
        
        ix = int(x)
        iy = int(y)
        iz = int(z)
        
        # Casting as int truncates. Must account for negative numbers.
        if ix <= x:
            x0 = ix
            x1 = ix + 1
        else:
            x0 = ix + 1
            x1 = ix
        
        if iy <= y:
            y0 = iy
            y1 = iy + 1
        else:
            y0 = iy + 1
            y1 = iy
        
        if iz <= z:
            z0 = iz
            z1 = iz + 1
        else:
            z0 = iz + 1
            z1 = iz
        
        return x0, y0, z0, x1, y1, z1
    
    def randomGradient(self, ix:int, iy:int, iz:int) -> tuple:
        """Return random gradient vector as a three-element tuple at a given
        Perlin lattice coordinate

        Uses spherical to rectangular coordinate conversion to generate vectors in 3D.
        Returned gradient vectors have a guaranteed length of sqrt(2) (i.e. rho = sqrt(2)). (this
        simplifies the normalization computation elsewhere in the code)
        

        Args:
            ix: integer x-coordinate of the Perlin lattice
            iy: integer y-coordinate of the Perlin lattice
            iz: integer z-coordinate of the Perlin lattice
        
        Returns:
            Pseudo-random gradient vector as a three-element tuple

        Notes:
            For other methods, see
            https://en.wikipedia.org/wiki/Perlin_noise*Implementation
        """
        
        theta = self._random_number_generator(ix, iy, iz, self._seed)/0xFFFFFFFF * 2*np.pi
        phi = self._random_number_generator(ix, iy, iz, self._seed+10)/0xFFFFFFFF * np.pi
        return (np.sqrt(2)*np.cos(theta)*np.sin(phi), np.sqrt(2)*np.sin(theta)*np.sin(phi), np.sqrt(2)*np.cos(phi))
    
    def _get_lattice_gradients(self, x:float, y:float, z:float) -> float:
        """Returns the 8 gradient vectors for the lattice unit containing the given point
            a,b,c,d,e,f,g,h refer to the points of the unit cube defined in _get_lattice_indices
        """
        x0, y0, z0, x1, y1, z1 = self._get_lattice_indices(x, y, z)
        va = self.randomGradient(x0, y0, z0)
        vb = self.randomGradient(x0, y1, z0)
        vc = self.randomGradient(x0, y0, z1)
        vd = self.randomGradient(x0, y1, z1)
        ve = self.randomGradient(x1, y1, z1)
        vf = self.randomGradient(x1, y0, z1)
        vg = self.randomGradient(x1, y1, z0)
        vh = self.randomGradient(x1, y0, z0)
        
        return va, vb, vc, vd, ve, vf, vg, vh
    
    def _lerp(self, t:float, a0:float, a1:float) -> float:
        """Linear interpolation between two values a0 and a1"""
        return (1-t)*a0 + t*a1

    def _cosine_interpolation(self, mu:float, y1:float, y2:float):
        mu2 = (1-np.cos(mu*np.pi))/2
        return (y1*(1-mu2)+y2*mu2)
    
    def _fade(self, t:float) -> float:
        """Fade function defined by Ken Perlin. Smooths out the linear interpolation by adjusting the t-value"""
        return ((6 * t - 15) * t + 10) * t**3
    
    def _fade_derivatives(self, t:float):
        """Returns the derivative of the above fade function. """
        return 30 * t * t * (t - 1) * (t - 1)
    
    def f(self, x:float, y:float, z:float) -> float:
        """Perlin noise function generator

        Args: 
            x: x-coordinate
            y: y-coordinate
            z: z-coordinate
        
        Returns: The value of the Perlin Noise function at the given coordinates
        """
        # x = round(x, self.precision)
        # y = round(y, self.precision)
        # z = round(z, self.precision)
        # if np.array([x, y, z]).tobytes() in self.perlin_values.keys():
        #     return self.perlin_values[np.array([x, y, z]).tobytes()]
        x0, y0, z0, x1, y1, z1 = self._get_lattice_indices(x,y,z)
        va, vb, vc, vd, ve, vf, vg, vh = self._get_lattice_gradients(x, y, z)
        
        # Vectors from point (x, y, z) to lattice points
        # a,b,c, etc. are same as in _get_lattice_indices
        # a = x0, y0, z0      e = x1, y1, z1
        # b = x0, y1, z0      f = x1, y0, z1
        # c = x0, y0, z1      g = x1, y1, z0
        # d = x0, y1, z1      h = x1, y0, z0
         
        ra = np.array([x - x0, y - y0, z - z0])
        rb = np.array([x - x0, y - y1, z - z0])
        rc = np.array([x - x0, y - y0, z - z1])
        rd = np.array([x - x0, y - y1, z - z1])
        re = np.array([x - x1, y - y1, z - z1])
        rf = np.array([x - x1, y - y0, z - z1])
        rg = np.array([x - x1, y - y1, z - z0])
        rh = np.array([x - x1, y - y0, z - z0])
        
        # Dot product between gradients and lattice vectors
        da = np.dot(ra, va)
        db = np.dot(rb, vb)
        dc = np.dot(rc, vc)
        dd = np.dot(rd, vd)
        de = np.dot(re, ve)
        df = np.dot(rf, vf)
        dg = np.dot(rg, vg)
        dh = np.dot(rh, vh)        

        # Adjusted interpolation weights based on the fade function
        u = self._fade((x-x0)/abs(x1-x0))
        v = self._fade((y-y0)/abs(y1-y0))
        w = self._fade((z-z0)/abs(z1-z0))
        
        # Linear interpolation 
        out = self._lerp(w, 
                         self._lerp(v, self._lerp(u, da, dh), self._lerp(u, db, dg)),
                         self._lerp(v, self._lerp(u, dc, df), self._lerp(u, dd, de)))
        # out = self._cosine_interpolation(w, 
        #                  self._cosine_interpolation(v, self._cosine_interpolation(u, da, dh), self._cosine_interpolation(u, db, dg)),
        #                  self._cosine_interpolation(v, self._cosine_interpolation(u, dc, df), self._cosine_interpolation(u, dd, de)))
        self.values.append(out)
        return out