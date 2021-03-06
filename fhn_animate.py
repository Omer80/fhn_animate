# -*- coding: utf-8 -*-
"""
#  Created on Sun Dec 27 11:52:01 2015
#
#  fhn_rhs.py
#  
#  Copyright 2015 Omer Tzuk <omertz@post.bgu.ac.il>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  FHN model
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from fhn_rhs import rhs
import scipy.integrate as integrate
import argparse


class FHNmodel:
    def __init__(self, e,a1,a0,d,l,n):
        self.l = float(l)
        self.n = int(n)
        self.dx = float(l/n)
        self.dx2 = self.dx**2
        self.parameters = {'e':e,'a1':a1,'a0':a0,'d':d}
        self.time_elapsed = 0
        self.state = self.init_state(n)
    def init_state(self,n):
        self.u0 = np.random.random((n,n))
        self.v0 = np.random.random((n,n))
        return np.ravel([self.u0,self.v0])
    def dstate_dt(self,y,t):
        vec_u , vec_v = np.split(y,2)
        u,v=vec_u.reshape(self.n,self.n),vec_v.reshape(self.n,self.n)
        u_t,v_t = rhs(u,v,self.parameters,self.dx)
        return np.ravel([u_t,v_t])
    def step(self, dt):
        """execute one time step of length dt and update state"""
        self.state = integrate.odeint(self.dstate_dt, self.state, [0, dt])[1]
        self.time_elapsed += dt
        
def main(args):
    fig = plt.figure()
    #------------------------------------------------------------
    # set up initial state and global variables
    fhn = FHNmodel(1,2,3,4,3,3)
    dt = 1./30 # 30 fps
    duration = 10
    #------------------------------------------------------------        
    
    # each frame
    ims = []
    im = plt.imshow(fhn.u0, cmap='viridis', animated=True)
    ims.append([im])
    time = 0
    while time < duration:
        fhn.step(dt)
        vec_u , vec_v = np.split(fhn.state,2)
        u = vec_u.reshape(fhn.n,fhn.n)
        im = plt.imshow(u, cmap='viridis', animated=True)
        ims.append([im])
        time+=dt
    
    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                    repeat_delay=1000)
    if args.filename is not None:
        print "Saving file.."
        ani.save('%s.mp4'%args.filename)
    ani.save('%s.mp4'%args.filename)
#    plt.show()
    
def add_parser_arguments(parser):
    "Add command-line arguments."""
    parser.add_argument('filename', type=str, nargs='?', default=None,
                        help='output movie file')
    return parser
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    args = parser.parse_args()
    main(args)