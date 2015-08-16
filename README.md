
To run an example, do the following::

1) Download surface data and domain data
mkdir clm-netcdf
svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/share/domains/domain.clm/domain.lnd.fv1.9x2.5_USGS.110713.nc ./clm-netcdf
svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/surfdata_map/surfdata_1.9x2.5_simyr2000_c141219.nc ./clm-netcdf

2) Run the example provied in 82x1_sparse_grid in MATLAB
>>CLM45SparseGridDriver('82x1_sparse_grid/82x1_sparse_grid.cfg')

