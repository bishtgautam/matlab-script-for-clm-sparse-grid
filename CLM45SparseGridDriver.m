% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% - Creates a CLM45 surface dataset and domain netcdf files in an
%   unstrustured grid format for a list sites given by latitude/longitude.
% 
% - The script uses already existing CLM45 surface datasets and create
%   new dataset by finding nearest neighbor for each site.
%
% INPUT:
%       fname = Configuration file name.
%
% EXAMPLE:
%  # Download data
%  mkdir clm-netcdf
%  svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/surfdata_map/surfdata_1.9x2.5_simyr2000_c141219.nc clm-netcdf/surfdata_1.9x2.5_simyr2000_c141219.nc
%  svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/share/domains/domain.clm/domain.lnd.fv1.9x2.5_USGS.110713.nc clm-netcdf/domain.lnd.fv1.9x2.5_USGS.110713.nc 
%
%  # Run matlab script
%  CLM45SparseGridDriver('82x1_sparse_grid/82x1_sparse_grid.cfg')
%
% Gautam Bisht (gbisht@lbl.gov)
% 05-28-2015
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [fsurdat, fdomain] = CLM45SparseGridDriver(cfg_filename,mesh)

disp('1) Reading configuration file')
cfg = ReadConfigurationFile(cfg_filename);

disp('2) Reading latitude/longitude @ cell centroid')
[lat,lon]   = ReadLatLon(cfg.site_latlon_filename);

disp('3) Computing latitude/longitude @ cell vertex')
[latv,lonv] = ComputeLatLonAtVertex(lat,lon, cfg.dlat, cfg.dlon);

disp('4) Creating CLM surface dataset')
fsurdat = CreateCLMUgridSurfdatForCLM45(lat, lon, ...
                              cfg.clm_gridded_surfdata_filename, ...
                              cfg.out_netcdf_dir, ...
                              cfg.clm_usrdat_name, ...
                              cfg.set_natural_veg_frac_to_one);

disp('5) Creating CLM domain')
%mesh = create_mpas_mesh_information();
fdomain = CreateCLMUgridDomainForCLM45(lat, lon, ...
                             latv, lonv, ...
                             cfg.clm_gridded_domain_filename, ...
                             cfg.out_netcdf_dir, ...
                             cfg.clm_usrdat_name, ...
                             mesh);

if (~isempty(cfg.landuse_timeseries_filename))
    disp('6) Creating CLM landuse timeseries file')
    flanduse = CreateCLMUgridLanduseTimeseries(lat, lon, ...
                                 cfg.landuse_timeseries_filename, ...
                                 cfg.out_netcdf_dir, ...
                                 cfg.clm_usrdat_name, ...
                                 cfg.set_natural_veg_frac_to_one);
end
