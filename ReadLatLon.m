% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Read latitude/longitude for site
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [lat,lon] = ReadLatLon(fname)

tmp_str = strsplit(fname,'.');

switch tmp_str{end}
    case 'txt'
        [lat,lon] = ReadLatLonFromTxt(fname);
    otherwise
        error('Unsupported format to read site level lat/lon');
end

end

function [lat,lon] = ReadLatLonFromTxt(fname)

fid = fopen(fname,'r');
if (fid == -1)
    error(['Unable to open file: ' fname])
end

npts = fscanf(fid,'%d',1);

lon = zeros(npts,1);
lat = zeros(npts,1);

for ii = 1:npts
    lat(ii,1) = fscanf(fid,'%f',1);
    lon(ii,1) = fscanf(fid,'%f',1);
end

fclose(fid);

end
