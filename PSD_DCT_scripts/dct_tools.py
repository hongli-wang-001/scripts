######################################################################
#
# Functions to perform fast discrete cosine and sine transforms and
# their inverses in one and two dimensions.  These functions work by
# wrapping the DFT function from numpy, rather than explicitly
# performing the cosine and sine transforms themselves.  The sine
# transforms take arrays whose first element is zero and return arrays
# whose first element is also zero.  This differs from some other
# implementations, which drop the first element, since it is always
# zero.
#
#   dct(y): Type-II discrete cosine transform (DCT) of real data y
#   idct(a): Type-II inverse DCT of a
#   dct2(y): 2D DCT of 2D real array y
#   idct2(a): 2D inverse DCT real array a
#
# Written by Mark Newman <mejn@umich.edu>, June 24, 2011
# You may use, share, or modify this file freely
#
######################################################################


from numpy import empty,arange,exp,real,imag,pi
from numpy.fft import rfft,irfft
import numpy as np

######################################################################
# 1D DCT Type-II

def dct(y):
    N = len(y)
    y2 = empty(2*N,float)
    y2[:N] = y[:]
    y2[N:] = y[::-1]

    c = rfft(y2)
    phi = exp(-1j*pi*arange(N)/(2*N))
    return real(phi*c[:N])


######################################################################
# 1D inverse DCT Type-II

def idct(a):
    N = len(a)
    c = empty(N+1,complex)

    phi = exp(1j*pi*arange(N)/(2*N))
    c[:N] = phi*a
    c[N] = 0.0
    return irfft(c)[:N]


######################################################################
# 2D DCT

def dct2(y):
    M = y.shape[0]
    N = y.shape[1]
    a = empty([M,N],float)
    b = empty([M,N],float)

    for i in range(M):
        a[i,:] = dct(y[i,:])
    for j in range(N):
        b[:,j] = dct(a[:,j])

    return b


######################################################################
# 2D inverse DCT

def idct2(b):
    M = b.shape[0]
    N = b.shape[1]
    a = empty([M,N],float)
    y = empty([M,N],float)

    for i in range(M):
        a[i,:] = idct(b[i,:])
    for j in range(N):
        y[:,j] = idct(a[:,j])

    return y


#######################################################################
# 2D Data Detrend

def detrend2d(DATA):

	# Detrend removes the trend from the data
	# Adapted by Sam Degelia from Matlab source
	nx = np.shape(DATA)[0]
	ny = np.shape(DATA)[1]
	x, y = np.meshgrid(range(1,ny+1), range(1,nx+1))
	Xcolv = np.reshape(x, (nx*ny,1))
	Ycolv = np.reshape(y, (nx*ny,1))
	Zcolv = np.reshape(DATA, (nx*ny,1))
	Const = np.empty_like(Zcolv)
	Const[:] = 1
	Coefficients = np.linalg.lstsq(np.concatenate((Xcolv, Ycolv, Const),axis=1), Zcolv)[0]
	Xcoeff = Coefficients[0]
	Ycoeff = Coefficients[1]
	Ccoeff = Coefficients[2]
	DATA_p = Xcoeff * x + Ycoeff * y + Ccoeff
	DATA_f = DATA - DATA_p
	return DATA_f

