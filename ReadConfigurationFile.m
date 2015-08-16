% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Reads a configuration file.
%
% INPUT:
%       fname = Configuration file name.
%
% Gautam Bisht (gbisht@lbl.gov)
% 05-28-2015
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function cfg = ReadConfigurationFile(fname)

% Initialization
cfg.site_latlon_filename           = '';
cfg.clm_gridded_surfdata_filename  = '';
cfg.clm_gridded_domain_filename    = '';
cfg.clm_usrdat_name                = '';
cfg.dlat                           = 0;
cfg.dlon                           = 0;
cfg.lon_min                        = -999;
cfg.lon_max                        = -999;

% Read the file
fid = fopen (fname,'r');

if (fid < 0)
    error(['Unable to open file: ' fname]);
end

while ~feof(fid)
    line = fgetl(fid);
    
    if (~isempty(line))
        if (strcmp(line(1),'%') == 0)
            tmp_string = strsplit(line,' ');
            if (length(tmp_string) ~= 2)
                error(['Incorrect enrty: ' line])                
            end
            
            switch lower(tmp_string{1})
                case 'site_latlon_filename'
                    cfg.site_latlon_filename = tmp_string{2};
                case 'clm_gridded_surfdata_filename'
                    cfg.clm_gridded_surfdata_filename = tmp_string{2};
                case 'clm_gridded_domain_filename'
                    cfg.clm_gridded_domain_filename = tmp_string{2};
                case 'clm_usrdat_name'
                    cfg.clm_usrdat_name = tmp_string{2};
                case 'dlat'
                    cfg.dlat = str2double(tmp_string{2});
                case 'dlon'
                    cfg.dlon = str2double(tmp_string{2});
                case 'lon_min'
                    cfg.lon_min = str2double(tmp_string{2});                    
                case 'lon_max'
                    cfg.lon_max = str2double(tmp_string{2});
                otherwise
                    error(['Unknown variable: ' tmp_string{1}])
            end
        end
    end
    
end

fclose(fid);

% Do some error checking
if (isempty(cfg.site_latlon_filename))
    error(['Stopping because entry for site_latlon_filename not found in ' fname])
end

if (isempty(cfg.clm_gridded_surfdata_filename))
    error(['Stopping because entry for clm_gridded_surfdata_filename not found in ' fname])
end

if (isempty(cfg.clm_gridded_domain_filename))
    error(['Stopping because entry for clm_gridded_domain_filename not found in ' fname])
end

if (isempty(cfg.clm_usrdat_name))
    error(['Stopping because entry for clm_usrdat_name not found in ' fname])
end

if (cfg.dlat == 0)
    error(['Stopping because entry for dlat not found in ' fname])
end

if (cfg.dlon == 0)
    error(['Stopping because entry for dlon not found in ' fname])
end

if (cfg.lon_min == -999)
    error(['Stopping because entry for lon_min not found in ' fname])
end

if (cfg.lon_max == -999)
    error(['Stopping because entry for lon_max not found in ' fname])
end

if (cfg.dlat < 0)
    error(['Stopping because dlat is negative in ' fname])
end

if (cfg.dlon < 0)
    error(['Stopping because dlon is negative in ' fname])
end

if ~(cfg.lon_min == 0 && cfg.lon_max == 360)
    error(['Stopping because lon_min and lon_max not equal to 0 and 360.'])
end


loc = find(cfg.site_latlon_filename == '/');
if (~isempty(loc))
    cfg.out_netcdf_dir = cfg.site_latlon_filename(1:loc(end)-1);
else
    cfg.out_netcdf_dir = './';
end



