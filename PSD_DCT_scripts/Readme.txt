1. there are two python module packages of Discrete Cosine Transform (DCT) -- dcst(including dct.py and dcts.pyc) and dct_tools(dct_tools.py and dct_tools.pyc). They are both developed by same author, and similar. User may import either one of the two modules in your script.

2. rtma3d_anl_2019111112.nc is analysis of RTMA3D for case of 2019-11-11_12:00UTC (netcdf format) for testing your code. To save space, there are only 4 2-D variables (TH2/Q2/U10/V10), and XLAT/XLONG for 2D plot.
   rtma3d_fgs_2019111112.nc is firstguess of RTMA3D for case of 2019-11-11_12:00UTC (netcdf format) for testing your code. To save space, there are only 4 2-D variables (TH2/Q2/U10/V10), and XLAT/XLONG for 2D plot.

3. psd_dct_2D_demo.py:
   1) in this python code, there is a function -- dct_psd_2d(u, nx_in, ny_in, dx_in, dy_in), which is the main function to apply the DCT to 2-D data array of "u".
     input arguments:
       2D array u is the 2-D variable in physical space, first dimension is Y, 2nd one is X (which is same as the 2-D array used in WRF-ARW model). nx_in is the dimension size (grid point number) in X-axix, ny_in is dimension size in Y-axis. dx_in is the grid space in X-direction in unit of meter, and dy_in is grid space in y-direction (e.g., for HRRR CONUS, dx_in=dy_in=3000.0 meters.).
     output arguments:
       total_power: 1-D array of power spectral density (unit depends on physical variable)
       k: 1-D array of wavenumbers (unit: radians/meter)
   2) the  rest part of this code is the main driver which reads in a 2-D field from a netcdf format file, then calls function dct_psd_2d to calculate the power spectral density (PSD) of the 2-D variable, and finally makes the plot of PSD.

   3) The algorithm in function dct_psd_2d to bin the 2D wavenumbers and convert to 1D wavenumberis not optimial, will be modified in later version.

4. scaledecomp_dct_2D_demo.py:
   python code to implement the scale seperation of the 2D data field (2D analysis increment in this demo code). This code is based on python code from Dr.Aaron Johnson(OU/SoM, OU/MAPS) with minor modifications. The contour plot was changed to color mesh plot.

5. Question, suggestion, comment are welcome. -->Gang Zhao, IMSG at NOAA/NWS/NCEP/EMC, email--> Gang.Zhao@noaa.gov
