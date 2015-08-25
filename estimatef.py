#!/usr/bin/python

#
# MODULE:    estimatef.py
# PURPOSE:   Estimating the fraction of sites under selection
#  


import sys
import numpy as np
import scipy.interpolate as interp
import operator as op


# input SFS
# sfs is 1-D arrary of site frequency(1,2,..., n-1)
#

def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom


def theta_w(sfs):
	N = len( sfs ) + 1
	a1 = np.array(xrange(1, N ))  
	return sum( sfs)/sum(1.0/a1) # sum 1 ~ N-1 class


def theta_pi(sfs):
	N = len( sfs )  + 1
	s = sum( i*( N - i )*sfs[ i -1 ] for i in range( 1, N) )
	return s/float( ncr(N, 2 ) ) 



def diff_f( test, ref ):
	f = 0.0
	# optional: smoothing
	n = len( test )
	x = np.array(range(1, len(  test)  +1))
	best_s = int( n * np.var(test))
	i_test = interp.UnivariateSpline (x, test, s = best_s)
	best_s = int( n * np.var(ref))
	i_ref = interp.UnivariateSpline (x, ref , s= best_s)
	sfs_test = i_test( x )
	sfs_ref  = i_ref( x )
	# scaling SFS_test by ratio of theta_pi
	a1 = theta_pi( sfs_test)/theta_pi( sfs_ref )
	selectiontype = "negative "
	f = ( 1 - theta_w(sfs_ref*a1)/theta_w(sfs_test) )
	if  f < 0.:
		# scale by theta_w
		selectiontype = "positive"
		s = np.asarray( [i*(n+1 - i) for i in range(1, n + 1 )] )/float(ncr(n+1, 2 ))
		a2 = theta_w( sfs_test)/theta_w(sfs_ref )
		f = -( 1 - theta_pi( ref*a2) /theta_pi( test) )
	return selectiontype , f 


#
# load SFS from files into array
#
if len( sys.argv ) == 3 :
	file_testsfs = open( sys.argv[1] )
	file_refsfs = open( sys.argv[2])
	testsfs, refsfs = (), ()
	for l in file_testsfs:
		testsfs = np.array([ int(x) for x in l.split()] )
	for l in file_refsfs:
		refsfs = np.array([ int(x) for x in l.split() ]) 
	print  "Fraction of sites under {} selection: {}".format( *diff_f( testsfs, refsfs )[0:2] )



