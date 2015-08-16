% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Computes latitude/longitude at vertices of grid cells.
% Only supports quadrilateral grid cells.
%
% INPUT:
%       lon  = Vector containing longitude @ cell-center.
%       lat  = Vector containing latitude @ cell-center.
%       dlat = Latitudinal grid spacing
%       dlon = Longitudinal grid spacing
%
% OUTPUT:
%       latv = Vector containing latitude @ cell-vertex.
%       lonv = Vector containing longitude @ cell-vertex.
%
% Gautam Bisht (gbisht@lbl.gov)
% 05-28-2015
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [latv,lonv] = ComputeLatLonAtVertex(lat, lon, dlat, dlon)

npts = length(lat);

latv = zeros(npts,4);
lonv = zeros(npts,4);

for ii = 1:npts
    lonv(ii,:) = [lon(ii)-dlon/2 lon(ii)+dlon/2 lon(ii)+dlon/2 lon(ii)-dlon/2];
    latv(ii,:) = [lat(ii)-dlat/2 lat(ii)-dlat/2 lat(ii)+dlat/2 lat(ii)+dlat/2];
end

