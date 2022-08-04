% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Creates an unstructured domain netCDF file for CLM45.
%
% INPUT:
%       lon_region = Vector containing longitude @ cell-center.
%       lat_region = Vector containing latitude @ cell-center.
%       lonv_region = Vector containing longitude @ cell-vertex.
%       latv_region = Vector containing latitude @ cell-vertex.
%       clm_gridded_domain_filename = Default CLM domain file name
%       out_netcdf_dir = Directory where CLM surface dataset will be saved
%       clm_usrdat_name = User defined name for CLM dataset
%
% Gautam Bisht (gbisht@lbl.gov)
% 01-02-2014
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function fname_out = CreateCLMUgridDomainForCLM45(lat_region, lon_region, ...
                                      latv_region, lonv_region, ...
                                      clm_gridded_domain_filename, ...
                                      out_netcdf_dir, ...
                                      clm_usrdat_name, mesh)


fname_out = sprintf('%s/domain_%s_%s.nc',out_netcdf_dir,clm_usrdat_name,datestr(now, 'cyymmdd'));
disp(['  domain: ' fname_out])

% Check if the file is available
[s,~]=system(['ls ' clm_gridded_domain_filename]);

if (s ~= 0)
   error(['File not found: ' clm_gridded_domain_filename]);
end

ncid_inp = netcdf.open(clm_gridded_domain_filename,'NC_NOWRITE');
ncid_out = netcdf.create(fname_out,'NC_CLOBBER');

info_inp = ncinfo(clm_gridded_domain_filename);

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_inp);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%                           Define dimensions
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
nv = size(lonv_region,2);
ni = size(lonv_region,1);
nj = 1;

for ii = 1:ndims
    [dimname, ndim] = netcdf.inqDim(ncid_inp,ii-1);
    switch dimname
        case 'ni'
            ndim = ni;
        case 'nj'
            ndim = nj;
        case 'nv'
            ndim = nv;
    end
    dimid(ii) = netcdf.defDim(ncid_out,dimname,ndim);
end

ndims = ndims + 1; dimid(ndims) = netcdf.defDim(ncid_out,'nCells',mesh.nCells); dim_nCells = ndims-1;
ndims = ndims + 1; dimid(ndims) = netcdf.defDim(ncid_out,'maxEdges',mesh.maxEdges); dim_maxEdges = ndims-1;
ndims = ndims + 1; dimid(ndims) = netcdf.defDim(ncid_out,'nEdges',mesh.nEdges); dim_nEdges = ndims-1;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%                           Define variables
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for ivar = 1:nvars
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_inp,ivar-1);
    varid(ivar) = netcdf.defVar(ncid_out,varname,xtype,dimids);
    varnames{ivar} = varname;
    
    for iatt = 1:natts
        attname = netcdf.inqAttName(ncid_inp,ivar-1,iatt-1);
        attvalue = netcdf.getAtt(ncid_inp,ivar-1,attname);
        
        netcdf.putAtt(ncid_out,ivar-1,attname,attvalue);
    end
    
end

nvars_orig = nvars;
nvars = nvars + 1;
netcdf.defVar(ncid_out,'cellsOnCell','NC_INT',[dim_maxEdges dim_nCells]);
netcdf.putAtt(ncid_out,nvars-1,'unit','unitless');
netcdf.putAtt(ncid_out,nvars-1,'long_name','List of cells that neighbor each cell.');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'cellsOnEdgeUp','NC_INT',[dim_nEdges]);
netcdf.putAtt(ncid_out,nvars-1,'unit','unitless');
netcdf.putAtt(ncid_out,nvars-1,'long_name','Cell indices that saddle a given edge.');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'cellsOnEdgeDown','NC_INT',[dim_nEdges]);
netcdf.putAtt(ncid_out,nvars-1,'unit','unitless');
netcdf.putAtt(ncid_out,nvars-1,'long_name','Cell indices that saddle a given edge.');

nvars = nvars + 1; 
netcdf.defVar(ncid_out,'edgesOnCell','NC_INT',[dim_maxEdges dim_nCells]);
netcdf.putAtt(ncid_out,nvars-1,'unit','unitless');
netcdf.putAtt(ncid_out,nvars-1,'long_name','List of edges that border each cell.');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'nEdgesOnCell','NC_INT',[dim_nCells]);
netcdf.putAtt(ncid_out,nvars-1,'unit','unitless');
netcdf.putAtt(ncid_out,nvars-1,'long_name','Number of edges that border each cell.');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'dcEdge','NC_DOUBLE',dim_nEdges);
netcdf.putAtt(ncid_out,nvars-1,'unit','m');
netcdf.putAtt(ncid_out,nvars-1,'long_name','Length of each edge, computed as the distance between cellsOnEdge.');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'dvEdge','NC_DOUBLE',dim_nEdges);
netcdf.putAtt(ncid_out,nvars-1,'unit','m');
netcdf.putAtt(ncid_out,nvars-1,'long_name','Length of each edge, computed as the distance between verticesOnEdge.');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'areaCell','NC_DOUBLE',dim_nCells);
netcdf.putAtt(ncid_out,nvars-1,'unit','m^2');
netcdf.putAtt(ncid_out,nvars-1,'long_name','Area of each cell in the primary grid.');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'xCell','NC_DOUBLE',dim_nCells);
netcdf.putAtt(ncid_out,nvars-1,'unit','m');
netcdf.putAtt(ncid_out,nvars-1,'long_name','x axis position of all cell centers');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'yCell','NC_DOUBLE',dim_nCells);
netcdf.putAtt(ncid_out,nvars-1,'unit','m');
netcdf.putAtt(ncid_out,nvars-1,'long_name','y axis position of all cell centers');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'zCell','NC_DOUBLE',dim_nCells);
netcdf.putAtt(ncid_out,nvars-1,'unit','m');
netcdf.putAtt(ncid_out,nvars-1,'long_name','z axis position of all cell centers');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'cosCell','NC_DOUBLE',dim_nCells);
netcdf.putAtt(ncid_out,nvars-1,'unit','unitless');
netcdf.putAtt(ncid_out,nvars-1,'long_name','cosine of the vertical angle between cell plane and gravity vector');

nvars = nvars + 1;
netcdf.defVar(ncid_out,'cosEdge','NC_DOUBLE',dim_nEdges);
netcdf.putAtt(ncid_out,nvars-1,'unit','unitless');
netcdf.putAtt(ncid_out,nvars-1,'long_name','cosine of the angle between unit vector joining cell neighbors and unit normal to the edge');

varid = netcdf.getConstant('GLOBAL');
[~,user_name]=system('echo $USER');
netcdf.putAtt(ncid_out,varid,'Created_by' ,user_name(1:end-1));
netcdf.putAtt(ncid_out,varid,'Created_on' ,datestr(now,'ddd mmm dd HH:MM:SS yyyy '))

netcdf.endDef(ncid_out);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%                           Copy variables
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
nvars = nvars_orig;
for ivar = 1:nvars
    
    data = netcdf.getVar(ncid_inp,ivar-1);
    [varname,vartype,vardimids,varnatts] = netcdf.inqVar(ncid_inp,ivar-1);
    
    switch varname
        case 'xc'
            data = lon_region;
        case 'yc'
            data = lat_region;
        case 'xv'
            data = lonv_region;
        case 'yv'
            data = latv_region;
        case 'mask'
            data = ones(length(lon_region),1);
        case 'frac'
            data = ones(length(lon_region),1);
        case 'area'
            if (size(lonv_region,2) == 3)
                ax = lonv_region(:,1);
                ay = latv_region(:,1);
                bx = lonv_region(:,2);
                by = latv_region(:,2);
                cx = lonv_region(:,3);
                cy = latv_region(:,3);
                
                data = 0.5*(ax.*(by-cy) + bx.*(cy-ay) + cx.*(ay-by));
            elseif (size(lonv_region,2) == 4)
                data = (lonv_region(:,1) - lonv_region(:,2)) .* (latv_region(:,1) - latv_region(:,3));
            else
                error('Added area computation')
            end
    end
    netcdf.putVar(ncid_out,ivar-1,data);
end



netcdf.close(ncid_inp);
netcdf.close(ncid_out);

ncwrite(fname_out,'cellsOnCell',mesh.cellsOnCell');
ncwrite(fname_out,'cellsOnEdgeUp',mesh.cellsOnEdge(:,1)');
ncwrite(fname_out,'cellsOnEdgeDown',mesh.cellsOnEdge(:,2)');
ncwrite(fname_out,'edgesOnCell',mesh.edgesOnCell');
ncwrite(fname_out,'nEdgesOnCell',mesh.nEdgesOnCell');
ncwrite(fname_out,'dcEdge',mesh.dcEdge);
ncwrite(fname_out,'dvEdge',mesh.dvEdge);
ncwrite(fname_out,'areaCell',mesh.areaCell);

[nx,ny] = size(mesh.xc);
ncwrite(fname_out,'xCell',reshape(mesh.xc,nx*ny,1));
ncwrite(fname_out,'yCell',reshape(mesh.yc,nx*ny,1));
ncwrite(fname_out,'zCell',reshape(mesh.zc,nx*ny,1));
ncwrite(fname_out,'cosCell',reshape(mesh.vcos,nx*ny,1));
ncwrite(fname_out,'cosEdge',mesh.cosEdge);


