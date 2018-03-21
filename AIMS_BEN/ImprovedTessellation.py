#!/usr/bin/env python
# $Id: ImprovedTessellation.py
# Author: Daniel R. Reese <dreese@bison.ph.bham.ac.uk>
# Copyright (C) Daniel R. Reese and contributors
# Copyright license: GNU GPL v3.0
#
#   AIMS is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with AIMS.  If not, see <http://www.gnu.org/licenses/>.
#

"""
This implements a :py:class:`ImprovedTessellation`, a subclass of 
:py:class:`scipy.spatial.Delaunay` which provides support for extrapolation
outside the grid of points, as well as various utility methods.
"""

__docformat__ = 'restructuredtext'

import math
import numpy as np
from scipy.spatial import Delaunay

class ImprovedTessellation(Delaunay):
    """
    A class which extends scipy.spatial.Delaunay thereby providing an improved
    tessellation which allows the user to find the nearest simplex outside the
    convex hull of the grid of points upon which the tessellation is built. 
    This facilitates extrapolation outside the grid of points.
    """

    def __init__(self,points):
        """
        :param points: coordinates of the points to triangulate
        :type points: (npoints, ndim) numpy float array
        """

        # most of the hard work is done here:
        Delaunay.__init__(self,points)

        # obtain the number of dimensions
        ndim = self.points.shape[1]

        # obtain the number of facets in the convex hull:
        nfacets = self.convex_hull.shape[0]

        # obtain the number of simplices:
        try:
            nsimplices = self.simplices.shape[0]
        except AttributeError:
            # this is necessary for outdated versions of scipy which
            # use the incorrect term "vertices" instead of "simplices"
            self.simplices = self.vertices
            nsimplices = self.simplices.shape[0]

        # create array representation of each facet in the convex hull:
        self.convex_hull_array = np.empty((nfacets,ndim,ndim),dtype=self.points.dtype)
        """Array representation of the facets which make up the convex hull"""
        for i in xrange(nfacets):
            for j in xrange(ndim):
                self.convex_hull_array[i,j,:] = self.points[self.convex_hull[i,j]]

        # create lookup array for facets:
        self.facet_to_simplex = np.empty((nfacets),dtype=self.simplices.dtype)
        """
        Lookup array which gives the index of the simplex corresponding to each facet
        from the convex hull
        """
        for i in xrange(nfacets):
            for j in xrange(nsimplices):
                if (set.issubset(set(self.convex_hull[i]),set(self.simplices[j]))):
                    self.facet_to_simplex[i] = j
                    break

    def find_simplex(self,pt):
        """
        Find the simplices containing the given points.  If the given point is outside
        the convex hull of the grid of points, find nearest simplex with the most
        favourable height.

        :param pt: point for which we're searching for the containing or nearest simplex
        :type pt: array-like of floats

        :return: the index of the containing or nearest simplex
        :rtype: int
        """

        ndim = self.points.shape[1]
        pt2 = np.array(pt,dtype=self.points.dtype)
        pt2 = pt2.reshape((1,ndim))
        val = Delaunay.find_simplex(self,pt2)
        if (val[0] != -1): return val

        # If we've made it here, it means the point is outside the convex hull
        # of the grid of point.  We therefore need to find the closest facet,
        # with the most favourable height.

        # find nearest facets from the convex hull:
        # NOTE: there can be more than one
        pt2 = pt2.reshape((ndim,))
        nfacets = self.convex_hull.shape[0]
        dmin = np.inf
        facets = []
        for i in xrange(nfacets):
            d = ImprovedTessellation.point_simplex_distance(pt2,self.convex_hull_array[i])
            if (abs(d-dmin) < 1e-10):
                dmin = min(d,dmin)
                facets.append(i)
                continue

            if (d < dmin):
                dmin = d
                facets = [i]

        # choose the facet with most favourable height:
        if (len(facets) > 1):
            heightmax = 0.0
            simplex_array = np.empty((ndim+1,ndim),dtype=self.points.dtype)
            for i in facets:
                simplex_array[0:ndim,:] = self.convex_hull_array[i]
                simplex_array[ndim,:] = pt2
                # this gives the height divided by ndim:
                height = ImprovedTessellation.simplex_volume(simplex_array) \
                       / ImprovedTessellation.simplex_volume(self.convex_hull_array[i])
                if (height > heightmax):
                    heightmax = height
                    val = [self.facet_to_simplex[i]]
        else:
            val = [self.facet_to_simplex[facets[0]]]

        return val

    def find_barycentres(self):
        """
        Find the barycentres of the simplices.
        
        :return: the set of barycentres
        :rtype: 2D numpy float array
        """
    
        nsimplices = self.simplices.shape[0]
        ndim = self.points.shape[1]
        
        result = np.zeros((nsimplices,ndim),dtype=self.points.dtype)
        for i in xrange(nsimplices):
            for j in self.simplices[i,:]:
                result[i,:] += self.points[j,:]
        result *= (1.0/(ndim+1))
    
        return result

    def find_barycentre_distances(self):
        """
        Find minimun and maximum distances between the barycentres 
        of the simplices and points which make up these simplices.
    
        :return: the set of distances
        :rtype: 2D numpy float array
        """
    
        nsimplices = self.simplices.shape[0]
        ndim = self.points.shape[1]
    
        result = np.empty((nsimplices,2),dtype=self.points.dtype)
        result[:,0] = np.inf
        result[:,1] = 0.0
        barycentres = self.find_barycentres()
        for i in xrange(nsimplices):
            for j in self.simplices[i,:]:
                diff = barycentres[i,:] - self.points[j,:]
                dist = np.dot(diff,diff)
                if (result[i,0] > dist): result[i,0] = dist
                if (result[i,1] < dist): result[i,1] = dist
        return np.sqrt(result)

    def simplex_to_simplex_array(self, val):
        """
        Convert an index representing a simplex to an array representation.

        :param val: index corresponding to a simplex
        :type val: integer
 
        :return: an array representation of the simplex
        :rtype: (ndim+1,ndim) numpy float array
        """

        assert (val >= 0), "Negative index not allowed for a simplex"
        assert (val < self.simplices.shape[0]), "Index larger than the number of simplices"
 
        ndim = self.points.shape[1]
        result = np.zeros((ndim+1,ndim),dtype=self.points.dtype)
        for i in xrange(ndim+1):
             result[i,:] = self.points[self.simplices[val,i],:]
    
        return result

    def find_distances(self,points):
        """
        Find minimun and maximum distances between a given set
        of points, and the points of the closest simplices.

        :param points: coordinates of the points for which we would like to
                       find the maximum and minimum distances
        :type points: (npoints, ndim) numpy float array

        :return: the set of distances
        :rtype: 2D numpy float array
        """

        assert (self.points.shape[1] == points.shape[1]), \
           "Dimension mismatch between input and tessellation points"

        npoints = points.shape[0]
        ndim    = points.shape[1]

        result = np.empty((npoints,2),dtype=self.points.dtype)
        result[:,0] = np.inf
        result[:,1] = 0.0
        for i in xrange(npoints):
            val = self.find_simplex(points[i,:])
            for j in self.simplices[val[0],:]:
                diff = points[i,:] - self.points[j,:]
                dist = np.dot(diff,diff)
                if (result[i,0] > dist): result[i,0] = dist
                if (result[i,1] < dist): result[i,1] = dist
        return np.sqrt(result)

    def find_simplex_volumes(self):
        """
        Find volumes of simplices.
    
        :return: the set of volumes
        :rtype: 1D numpy float array
        """
    
        nsimplices = self.simplices.shape[0]
        ndim = self.points.shape[1]
    
        result = np.zeros((nsimplices,),dtype=self.points.dtype)
        simplex_array = np.empty((ndim+1,ndim),dtype=self.points.dtype)
        for i in xrange(nsimplices):
            for j in xrange(ndim+1):
                simplex_array[j,:] = self.points[self.simplices[i,j]]
            result[i] = ImprovedTessellation.simplex_volume(simplex_array)
        return result

    @staticmethod
    def point_simplex_distance(pt,simplex):
        """
        Recursive function which finds the distance between a simplex and a point.
        Reference: `Golubitsky, O., Mazalov, V., Watt, S. M. (2015) <http://www.researchgate.net/publication/267801141__An_Algorithm_to_Compute_the_Distance_from_a_Point_to_a_Simplex>`_

        :param pt: n-dimensional point
        :type pt: numpy 1D float array

        :param simplex: a (d-1)-dimensional simplex (d-1 <= n)
        :type simplex: (d,n) numpy array

        :return: the distance between pt and the closest point in simplex
        :rtype: float
        """

        (d,n) = simplex.shape

        # sanity check:
        if (d > n+1): return np.nan
        if (d == 0): return np.nan

        if (d == 1):
            diff = pt - simplex[0,:]
            return math.sqrt(np.dot(diff,diff))

        d-=1
        mat = np.empty((d,d),dtype=np.float64)
        v   = np.empty((d,),dtype=np.float64)
        for i in xrange(d):
            diff = simplex[i+1,:]-simplex[0,:]
            v[i] = np.dot(diff,pt-simplex[0,:])
            for j in xrange(i,d):
                mat[i,j] = np.dot(diff,simplex[j+1,:]-simplex[0,:])
                mat[j,i] = mat[i,j]

        alp = np.linalg.solve(mat,v)

        for i in xrange(d):
            if (alp[i] < 0.0):
                simplex_reduced = np.empty((d,n),dtype=np.float64)
                for j in xrange(i+1): simplex_reduced[j,:] = simplex[j,:]
                for j in xrange(i+1,d): simplex_reduced[j,:] = simplex[j+1,:]
                return ImprovedTessellation.point_simplex_distance(pt,simplex_reduced)

        my_sum = np.sum(alp)
        if (my_sum > 1.0):
            simplex_reduced = np.empty((d,n),dtype=np.float64)
            for j in xrange(d): simplex_reduced[j,:] = simplex[j+1,:]
            return ImprovedTessellation.point_simplex_distance(pt,simplex_reduced)

        pt_bis = np.zeros((n,),dtype=np.float64)
        for i in xrange(d):
            pt_bis += alp[i]*simplex[i+1,:]
        pt_bis += (1.0-my_sum)*simplex[0,:]
        diff = pt_bis - pt
        return math.sqrt(np.dot(diff,diff))

    @staticmethod
    def simplex_volume(simplex):
        """
        Calculate the volume of a simplex. The approach used here also works
        for simplices embedded in higher dimensional spaces (e.g. a triangle
        in 3D space).

        :param simplex: a (d-1)-dimensional simplex (d-1 <= n)
        :type simplex: (d,n) numpy array

        :return: the volume of the simplex
        :rtype: float
        """

        (d,n) = simplex.shape

        # sanity check
        if (d > n+1): return np.nan
        if (d == 0): return np.nan

        # easy exit
        if (d == 1): return 1.0  # this result is based on induction

        # construct basis for simplex
        basis = np.empty((d-1,n),dtype=np.float64)
        for i in xrange(1,d):
            basis[i-1,:] = simplex[i,:]-simplex[0,:]

        # orthogonalise basis and get squared volume:
        result = 1.0
        norm2 = np.empty((d-1,),dtype=np.float64)
        for i in xrange(d-1):
            for j in xrange(i):
                coef = np.dot(basis[i],basis[j])/norm2[j]
                basis[i] -= basis[j]*coef
            norm2[i] = np.dot(basis[i],basis[i])
            result *= norm2[i]

        # find volume by taking square root and dividing by (d-1)!
        # (see wikipedia for an explanation of the 1/(d-1)! coefficient)
        return math.sqrt(result)/math.factorial(d-1)

