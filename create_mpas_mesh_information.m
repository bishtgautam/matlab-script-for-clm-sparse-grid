function mesh = create_mpas_mesh_information(nx,ny,X,Y)

nx = 5;
ny = 4;

X = 1000; % [m]
Y = 1000; % [m]

dx = X/nx;
dy = Y/ny;

nCells = nx*ny;
maxEdges = 4;

cellsOnCell  = zeros(nCells,maxEdges);
edgesOnCell  = zeros(nCells,maxEdges);
edge_ids     = zeros(nCells,maxEdges);
nEdgesOnCell = ones(nCells,1)*4;
areaCell     = zeros(nCells,1);
xCell        = zeros(nCells,1);
yCell        = zeros(nCells,1);
zCell        = zeros(nCells,1);
dcEdge       = zeros(nCells*maxEdges,1);
dvEdge       = zeros(nCells*maxEdges,1);

cell_ids     = reshape([1:nCells],nx,ny);
vert_ids     = reshape([1:(nx+1)*(ny+1)],nx+1,ny+1);


nEdges = 0;

for ii = 1:nx+1
    for jj = 1:ny+1
        
        xVertex(vert_ids(ii,jj)) = ii*dx;
        yVertex(vert_ids(ii,jj)) = jj*dy;
        zVertex(vert_ids(ii,jj)) = (ii-1)*(jj-1);
    end
end
    

for ii = 1:nx
    for jj = 1:ny
        
        cell_id = cell_ids(ii,jj);

        areaCell(cell_id) = dx*dy;
        xCell(cell_id) = mean(xVertex([vert_ids(ii,jj) vert_ids(ii+1,jj) vert_ids(ii+1,jj+1) vert_ids(ii,jj+1)]));
        yCell(cell_id) = mean(yVertex([vert_ids(ii,jj) vert_ids(ii+1,jj) vert_ids(ii+1,jj+1) vert_ids(ii,jj+1)]));
        zCell(cell_id) = mean(zVertex([vert_ids(ii,jj) vert_ids(ii+1,jj) vert_ids(ii+1,jj+1) vert_ids(ii,jj+1)]));
        
        count = 0;
        
        if (ii > 1)
            nEdges = nEdges+1;
            cellsOnEdge(nEdges,1:2) = [cell_ids(ii-1,jj) cell_ids(ii,jj)];
            %dcEdge(nEdges,1) = dx;
            %dvEdge(nEdges,1) = dy;

            count = count + 1;
            cellsOnCell(cell_id, count) = cell_ids(ii-1,jj);
            %edgesOnCell(cell_id, count) = nEdges;
        end
        
        iedge = 1;
        edge_idx = (cell_id-1)*maxEdges + iedge;
        edgesOnCell(cell_id,iedge) = edge_idx;
        dcEdge(edge_idx) = dx;
        dvEdge(edge_idx) = dy;
        
        if (ii < nx)
            count = count + 1;
            cellsOnCell(cell_id,count) = cell_ids(ii+1,jj);
        end
        iedge = 2;
        edge_idx = (cell_id-1)*maxEdges + iedge;
        edgesOnCell(cell_id,iedge) = edge_idx;
        dcEdge(edge_idx) = dx;
        dvEdge(edge_idx) = dy;
        
        if (jj > 1)
            nEdges = nEdges+1;
            cellsOnEdge(nEdges,1:2) = [cell_ids(ii,jj-1) cell_ids(ii,jj)];
            %dcEdge(nEdges,1) = dy;
            %dvEdge(nEdges,1) = dx;

            count = count + 1;
            cellsOnCell(cell_id,count) = cell_ids(ii,jj-1);
        end
        iedge = 3;
        edge_idx = (cell_id-1)*maxEdges + iedge;
        edgesOnCell(cell_id,iedge) = edge_idx;
        dcEdge(edge_idx) = dy;
        dvEdge(edge_idx) = dx;
        
        if (jj < ny)
            count = count + 1;
            cellsOnCell(cell_id,count) = cell_ids(ii,jj+1);
            %edgesOnCell(cell_id, count) = nEdges;
        end
        
        iedge = 4;
        edge_idx = (cell_id-1)*maxEdges + iedge;
        edgesOnCell(cell_id,iedge) = edge_idx;
        dcEdge(edge_idx) = dy;
        dvEdge(edge_idx) = dx;
        
    end
end

% for iedge = 1:nEdges
%     icell1 = cellsOnEdge(iedge,1);
%     icell2 = cellsOnEdge(iedge,2);
%     
%     loc = find(cellsOnCell(icell1,:) == icell2); edgesOnCell(icell1,loc) = iedge;
%     loc = find(cellsOnCell(icell2,:) == icell1); edgesOnCell(icell2,loc) = iedge;
% end


%cellsOnCell = reshape(cellsOnCell',nCells*maxEdges,1);

mesh.cellsOnCell = cellsOnCell;
mesh.edgesOnCell = edgesOnCell;
mesh.nEdgesOnCell = nEdgesOnCell;
mesh.cellsOnEdge = cellsOnEdge;
mesh.dcEdge = dcEdge;
mesh.dvEdge = dvEdge;
mesh.areaCell = areaCell;
mesh.nCells = nCells;
mesh.maxEdges = maxEdges;
mesh.nEdges = nCells*maxEdges;
mesh.cell_ids = cell_ids;