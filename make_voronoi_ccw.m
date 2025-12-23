function [newV, newC, newC_mat, adjCV, adjCC] = make_voronoi_ccw(cx, cy)
% newV (N x 2 matrix, where N is the number of cells) gives the list of vertex positions
% newC is a cell with N entries, each entry gives the list of vertex neighbors for a cell. The numbers in each entry correspond to the index in newV
% newC_mat is just newC in matrix format
% adjCV (N x 2*N matrix), this is the Cell-Vertex adjacency matrix. e.g. the (i,j) entry is 1 if cell i and vertex j are neighbors, 0 if they are not
% adjCC (N x N matrix), this is the Cell-Cell adjacency matrix. e.g. the (i,j) entry is 1 if cell i and cell j are neighbors, 0 if they are not



% Author: Max Bi
% e-mail address: bdpmax@gmail.com

% Parameters for voronoi tesselation
% Maximum number of neighbors a cell can have, chose 20 to be safe
maxzz = 20;
% Parameters for voronoi tesselation
% distance tolerance for distinguishing two vertices
tol=1e-8;



npts = length(cx);
box_size = sqrt(npts);
cx_original = cx;
cy_original = cy;

cx = [[cx-1; cx; cx+1]; [cx-1; cx; cx+1]; [cx-1; cx; cx+1]];
cy = [[cy-1; cy-1; cy-1]; [cy; cy; cy]; [cy+1; cy+1; cy+1]];


lmin = 0.2;
ind_include = cx > - lmin & cx < (1+lmin) & cy > - lmin & cy < (1+lmin);
cx = cx(ind_include);
cy = cy(ind_include);

ind_original =  cx > 0.0 & cx < 1.0 & cy > 0.0 & cy < 1.0;

[V, oldC]=voronoin([cx,cy]);
% newC = oldC((npts*4+1):(npts*5));
newC = oldC(ind_original);
newV_ind = unique([newC{:}]);
newV =mod(V(newV_ind,:),1);
[newV, ~, ind] = consolidator(newV,[],[],tol);

while size(newV,1) ~= npts*2
    %     disp(size(newV,1));
    cx = cx_original + randn(npts,1).*sqrt(10*tol);
    cy = cy_original + randn(npts,1).*sqrt(10*tol);
    
    cx = [[cx-1; cx; cx+1]; [cx-1; cx; cx+1]; [cx-1; cx; cx+1]];
    cy = [[cy-1; cy-1; cy-1]; [cy; cy; cy]; [cy+1; cy+1; cy+1]];
    
    [V, oldC]=voronoin([cx,cy]);
    newC = oldC((npts*4+1):(npts*5));
    newV_ind = unique([newC{:}]);
    newV =mod(V(newV_ind,:),1);
    [newV, ~, ind] = consolidator(newV,[],[],tol);
    
end

for i = 1:npts
    [~,loc]=ismember(newC{i},newV_ind);
    newC{i} = ind(loc)';
end




% Matrix version of newC
newC_mat = cell2mat(...
    cellfun(...
    @(xx)cat(2,xx,zeros(1,maxzz-length(xx))),...
    newC,'UniformOutput',false...
    )...
    );
% Cell-cell adjacency matrix
adjCC = zeros(npts,npts);
for i = 1:length(newV)
    [m, ~]=find(newC_mat==i);
    adjCC(m(1),m(2))=1;
    adjCC(m(2),m(1))=1;
    adjCC(m(1),m(3))=1;
    adjCC(m(3),m(1))=1;
    adjCC(m(2),m(3))=1;
    adjCC(m(3),m(2))=1;
end
% Ordered Cell-Vertex adjacency matrix
adjCV = zeros(npts,length(newV));
for i = 1:npts
    zz = length(newC{i});
    adjCV(i,newC{i}) = 1:zz;
end
newV = newV*sqrt(npts);

%
for i = 1:npts
    vlist = newC{i};
    zz = length(newC{i});
    
    Vx = newV(vlist,1);
    Vy = newV(vlist,2);
    % make sure all vertices are within the same PBC image
    Vx = Vx(end) + pbc_dist(Vx-Vx(end),box_size);
    Vy = Vy(end) + pbc_dist(Vy-Vy(end),box_size);
    %     [~,ord] = sort(mod(atan2(Vy-Vy(1),Vx-Vx(1)),2*pi));
    [~,ord] = sort(mod(atan2(Vy-mean(Vy),Vx-mean(Vx)),2*pi));
    
    vlist = vlist(ord);
    newC{i} = vlist;
    
end

end

