# -*- coding: utf-8 -*-
"""
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
#   FHN model
"""
__version__=1.0
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""
#from numbapro import jit, int32, float32
from numbapro import autojit
from derivatives import neumann_laplacian_2D_2nd_order as laplacian
@autojit(target="cpu")
def rhs(u,v,parameters,dx2):
    """ RHS of FHN model, only reaction terms """
    e = parameters['e']
    a1 = parameters['a1']
    a0 = parameters['a0']
    d = parameters['d']
    u_t = u - u**3 - v  + laplacian(u,dx2)
    v_t = e*(u - a1*v - a0) + d*laplacian(v,dx2)
    return u_t, v_t
