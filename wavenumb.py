from numpy import linspace , arange , reshape ,zeros
from scipy.fftpack import fft2 , fftfreq
from cmath import pi
import scipy.io as sc
import numpy as np
import idlsave

btot_rot=sc.readsav('magn_field/b_los_6302_aligned2ca.sav',verbose=True,python_dict=True)['btot_rot']
im=btot_rot[26,650:800,300:450]

# create some arbitrary data
#some_data = arange(0.0 , 16384.0 , dtype = complex)

some_data_grid = im
# reshape it to be a 128x128 2d grid
#some_data_grid = reshape(some_data , (128 , 128) )

# assign some real spatial co-ordinates to the grid points   
# first define the edge values
x_min = 0
x_max = im.shape[0]
y_min = 0
y_max = im.shape[0]

# then create some empty 2d arrays to hold the individual cell values
x_array = zeros( (im.shape[0],im.shape[0]) , dtype = float )
y_array = zeros( (im.shape[0],im.shape[0]) , dtype = float )

# now fill the arrays with the associated values
for row , y_value in enumerate(linspace (y_min , y_max , num = im.shape[0]) ):

  for column , x_value in enumerate(linspace (x_min , x_max , num = im.shape[0]) ):

    x_array[row][column] = x_value
    y_array[row][column] = y_value

# now for any row,column pair the x_array and y_array hold the spatial domain
# co-ordinates of the associated point in some_data_grid

# now use the fft to transform the data to the wavenumber domain
some_data_wavedomain = fft2(some_data_grid)

# now we can use fftfreq to give us a base for the wavenumber co-ords
# this returns [0.0 , 1.0 , 2.0 , ... , 62.0 , 63.0 , -64.0 , -63.0 , ... , -2.0 , -1.0 ]
n_value = fftfreq( im.shape[0] , (1.0 / im.shape[0] ) )

# now we can initialize some arrays to hold the wavenumber co-ordinates of each cell
kx_array = zeros( (im.shape[0],im.shape[0]) , dtype = float )
ky_array = zeros( (im.shape[0],im.shape[0]) , dtype = float )

# before we can calculate the wavenumbers we need to know the total length of the spatial
# domain data in x and y. This assumes that the spatial domain units are metres and
# will result in wavenumber domain units of radians / metre.
x_length = x_max - x_min
y_length = y_max - y_min

# now the loops to calculate the wavenumbers
for row in xrange(128):

  for column in xrange(128):

    kx_array[row][column] = ( 2.0 * pi * n_value[column] ) / x_length
    ky_array[row][column] = ( 2.0 * pi * n_value[row] ) / y_length
