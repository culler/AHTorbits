#   This file is part of the program Orbits.
#
#   Copyright (C) 2014-2017 by Marc Culler and others. 
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#   Project homepage: https://bitbucket.org/marc_culler/AHTorbits
#   Author homepage: http://marc-culler.info

from __future__ import print_function
from orbits import *

# Compute the number of components of a curve on a genus 2 surface
# split along a separating curve bounding two one-holed tori.
# We assume there are three parallels bands, of weights p, q, r, on each
# one-holed torus, and a shift of n on the circle.

def PG(p,q,r,n):
   N = p+q+r
   n = n%(2*N)
   G= Pseudogroup([ Flip([1,p],[p+q+r+1, 2*p+q+r]),
		    Flip([p+1, p+q], [2*p+q+r+1, 2*p+2*q+r]),
		    Flip([p+q+1, p+q+r], [2*p+2*q+r+1, 2*p+2*q+2*r]),
		    Flip([2*N+1,2*N+r],[2*N+p+q+r+1, 2*N+p+q+2*r]),
		    Flip([2*N+r+1, 2*N+r+q], [p+q+2*r+1, p+2*q+2*r]),
		    Flip([2*N+r+q+1, 2*N+p+q+r], [2*N+p+2*q+2*r+1, 2*N+2*p+2*q+2*r]),
		    Shift([1, 2*N - n],[2*N+n+1,4*N]),
                    Shift([2*N-n+1,2*N],[2*N+1, 2*N+n])], 
                Interval(1,4*N))
   return G.reduce()

if __name__ == '__main__':
   print('p = q = r = n = 1:')
   print(PG(1,1,1,1), 'components')
