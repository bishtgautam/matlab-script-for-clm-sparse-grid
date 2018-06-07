
## Objective

Scripts to create ELM/CLM surface dataset and domain NetCDF files for a list of user-specified grid cells.
The scripts create new datasets by extracting nearest neighbor data from existing global datasets.
The generated surface dataset and domain files are in the unstructured-grid format.

## Running the code

1. Download the code

```
git clone https://github.com/bishtgautam/matlab-script-for-clm-sparse-grid
```


2. Download surface data and domain data to run the 82x1_sparse_grid example

```
mkdir clm-netcdf
svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/share/domains/domain.clm/domain.lnd.fv1.9x2.5_USGS.110713.nc ./clm-netcdf
svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/surfdata_map/surfdata_1.9x2.5_simyr2000_c141219.nc ./clm-netcdf
```

3. Launch MATLAB

```matlab
<MATLAB_INSTALLATION_DIR>/bin/matlab -nodesktop
```

4. Run the example

```matlab
>> cd <matlab-script-for-clm-sparse-grid>
>> CLM45SparseGridDriver('82x1_sparse_grid/82x1_sparse_grid.cfg');
1) Reading configuration file
2) Reading latitude/longitude @ cell centroid
3) Computing latitude/longitude @ cell vertex
4) Creating CLM surface dataset
  surface_dataset: 82x1_sparse_grid/surfdata_82x1_sparse_grid_c180227.nc
  WARNING: Site with (lat,lon) = (43.730000,288.750000) has more than one cells that are equidistant.
           Picking the first closest grid cell.
		Possible grid cells: 44.526316 287.500000
		Possible grid cells: 44.526316 290.000000
5) Creating CLM domain
  domain: 82x1_sparse_grid/domain_82x1_sparse_grid_c180227.nc
```

## License

[BSD-3](https://github.com/bishtgautam/matlab-script-for-clm-sparse-grid/License.txt)
