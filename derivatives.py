import numpy as np
#from numbapro import autojit
"""
# Copyright (C) 2009 by 
#    Omer Tzuk    <omertz@post.bgu.ac.il>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY <copyright holder> ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
__version__=1.0
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""
#@autojit(target="cpu")
def neumann_laplacian_2D_2nd_order(u,dx2):
    """Return finite difference Laplacian approximation of 2d array.

    Uses Neumann boundary conditions and a 2nd order approximation.
    
    """
    laplacian = np.zeros(u.shape)
    laplacian[1:-1, 1:-1] =  (u[1:-1,2:] + u[1:-1,:-2]  
                              + u[2:,1:-1] + u[:-2,1:-1] 
                              -  4.0*u[1:-1,1:-1])/(dx2)
    # Neumann boundary conditions
    # edges
    laplacian[0,1:-1] =  (u[0,2:] + u[0,:-2]
                          + 2.0*u[1,1:-1] 
                          - 4.0*u[0,1:-1])/(dx2)
    laplacian[-1,1:-1] =  (u[-1,2:] + u[-1,:-2]
                           + 2.0*u[-2,1:-1] 
                           - 4.0*u[-1,1:-1])/(dx2)
    laplacian[1:-1,0] = (2.0*u[1:-1,1] 
                         + u[2:,0] + u[:-2,0] 
                         -  4.0*u[1:-1,0])/(dx2)
    laplacian[1:-1,-1] =   (2.0*u[1:-1,-2] 
                            + u[2:,-1] + u[:-2,-1] 
                            -  4.0*u[1:-1,-1])/(dx2)
    # corners
    laplacian[0,0]  =(2.0*u[0,1]   + 2.0*u[1,0]   - 4.0*u[0,0])/(dx2)
    laplacian[-1,0] =(2.0*u[-1,1]  + 2.0*u[-2,0]  - 4.0*u[-1,0])/(dx2)
    laplacian[0,-1] =(2.0*u[0,-2]  + 2.0*u[1,-1]  - 4.0*u[0,-1])/(dx2)
    laplacian[-1,-1]=(2.0*u[-1,-2] + 2.0*u[-2,-1] - 4.0*u[-1,-1])/(dx2)


    return laplacian

if __name__ == "__main__":
	pass
