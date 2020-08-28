function [ii_idx,jj_idx] = find_nearest_neighbor(lati_source,long_source,lati_dest,long_dest)

% find the index
for ii=1:size(long_dest,1)
    for jj=1:size(long_dest,2)
        dist = (long_source - long_dest(ii,jj)).^2 + (lati_source - lati_dest(ii,jj)).^2;
        [nearest_cell_i_idx, nearest_cell_j_idx] = find( dist == min(min(dist)));
        if (length(nearest_cell_i_idx) > 1)
            disp(['  WARNING: Site with (lat,lon) = (' sprintf('%f',lati_dest(ii,jj)) ...
                sprintf(',%f',long_dest(ii,jj)) ') has more than one cells ' ...
                'that are equidistant.' char(10) ...
                '           Picking the first closest grid cell.']);
            for kk = 1:length(nearest_cell_i_idx)
                disp(sprintf('\t\tPossible grid cells: %f %f', ...
                    lati_source(nearest_cell_i_idx(kk),nearest_cell_j_idx(kk)), ...
                    long_source(nearest_cell_i_idx(kk),nearest_cell_j_idx(kk))));
            end
        end
        ii_idx(ii,jj) = nearest_cell_i_idx(1);
        jj_idx(ii,jj) = nearest_cell_j_idx(1);
    end
end
