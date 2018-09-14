% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Performs checks on fractional cover.
%
% INPUT:
%       varname                = name of the netcdf variable
%       data                   = data values for the sparse grid
%       set_natural_veg_frac_to_one = Flag to override default fractional cover
%                                by setting PCT_NATVEG = 1 and all other
%                                fractional covers to 0
%
% Gautam Bisht (gbisht@lbl.gov)
% 02-27-2018
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function data = PerformFractionCoverCheck(varname, data, set_natural_veg_frac_to_one)

switch varname
    case {'PCT_URBAN'}        
        if (set_natural_veg_frac_to_one)
            data = data * 0;
        elseif (max(sum(data,2) > 0))
            disp(' ')
            disp(['Warning: ' varname ' is not 0. for all grids' ])
            disp('         If you wish to create surface dataset')
            disp('         with only natural vegetation, set')
            disp('         set_natural_veg_frac_to_one = 1 in CFG file')
            disp(' ')
        end
        
    case {'PCT_CROP','PCT_WETLAND','PCT_LAKE','PCT_GLACIER'}
        if (set_natural_veg_frac_to_one)
            data = data * 0;
        elseif (max(data) > 0)
            disp(' ')
            disp(['Warning: ' varname ' is not 0. for all grids' ])
            disp('         If you wish to create surface dataset')
            disp('         with only natural vegetation, set')
            disp('         set_natural_veg_frac_to_one = 1 in CFG file')
            disp(' ')
        end
    case {'PCT_NATVEG'}
        if (set_natural_veg_frac_to_one)
            data = data*0 + 100;
        elseif (min(data) < 100)
            disp(' ')
            disp(['Warning: ' varname ' is not 100. for all grids' ])
            disp('         If you wish to create surface dataset')
            disp('         with only natural vegetation, set')
            disp('         set_natural_veg_frac_to_one = 1 in CFG file')
            disp(' ')
        end
end
