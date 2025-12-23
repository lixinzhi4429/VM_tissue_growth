function [newV, newC] = make_cell_division_new(cell_id, oldV, oldC)

V_count = size(oldV,1);
Ni = size(oldC,1);
vlist = oldC{cell_id};
z = length(vlist);
vx = oldV(vlist,1);
vy = oldV(vlist,2);
cx = mean(vx);
cy = mean(vy);


dx = repmat(vx, 1, z) - repmat(vx.', z, 1);
dy = repmat(vy, 1, z) - repmat(vy.', z, 1);
dr = sqrt(dx.^2+dy.^2);

dr_tru = triu(dr);

[maxValue, ~] = max(dr_tru(:));
[row_id, col_id] = find(dr_tru == maxValue);

% slope of the main axis through cenn center
slope = dy(row_id,col_id)/dx(row_id,col_id);

vlist_cylic = [vlist,vlist(1)];
XY1 = [];
for j = 1:z
    XY1 = [XY1; oldV(vlist_cylic(j),1) oldV(vlist_cylic(j),2) oldV(vlist_cylic(j+1),1) oldV(vlist_cylic(j+1),2)];
end

% the line perperdicualt to the main axis
slope_perp = -1/slope;
b_perp = cy - cx*slope_perp;

xmin = min(oldV(vlist,1))-0.1;
xmax = max(oldV(vlist,1))+0.1;

y1 = slope_perp*xmin + b_perp;
y2 = slope_perp*xmax + b_perp;
XY2 = [xmin,y1,xmax,y2];

out = lineSegmentIntersect(XY1,XY2);

intX = out.intMatrixX;
intY = out.intMatrixY;
% ids = intX~=0; % in rare cases, not working if X is 0 only
ids = logical((intX~=0)+(intY~=0));
intX = intX(ids);
intY = intY(ids);


% this part not working well
% % the intersection points might coincide with existing vertices, which will
% % case NaN in force calcualtion, perturb them a little
% int_V1 = [intX(1),intY(1)];
% int_V2 = [intX(2),intY(2)];
% 
% dd1 = abs(int_V1 - oldV(vlist,1));
% dd2 = abs(int_V2 - oldV(vlist,1));
% 
% temp1 = find(dd1<1e-10);
% temp2 = find(dd2<1e-10);
% 
% if ~isempty(temp1) 
%     int_V1 = int_V1 + 0.01;
% elseif ~isempty(temp2)
%     int_V2 = int_V2 + 0.01;
% end
% 
% 
% intX = [int_V1(1);int_V2(1)];
% intY = [int_V1(2);int_V2(2)];


[~,~, adjCV, ~, ~, ~, ~] = make_adj_fbc(oldV, oldC);


inter_matrix = out.intAdjacencyMatrix;
inter_edges = find(inter_matrix);

v1 = vlist_cylic(inter_edges(1));
v2 = vlist_cylic(inter_edges(1)+1);
cneibs1 = find(adjCV(:,v1) & adjCV(:,v2));
cneib1 = cneibs1(cneibs1~=cell_id);


v3 = vlist_cylic(inter_edges(2));
v4 = vlist_cylic(inter_edges(2)+1);
cneibs2 = find(adjCV(:,v3) & adjCV(:,v4));
cneib2 = cneibs2(cneibs2~=cell_id);

vertices_split = [v1,v2,v3,v4];
% the two intersections add two new vertices

id1 = inter_edges(1);
id2 = inter_edges(2);

id_min = min([id1,id2]);
id_max = max([id1,id2]);
if id_max == z
    vlist_temp = [vlist(1:id_min),V_count+1,vlist(id_min+1:z),V_count+2];
else
    vlist_temp = [vlist(1:id_min),V_count+1,vlist(id_min+1:id_max),V_count+2,vlist(id_max+1:z)];
end

% make the vertices in CCW order
% errors occur for concave cell shapes, how to avoid?


% don't use vertex position to order, it fails for concave cell shapes.
% insert new vertex indices into old cell list correspondingly

% [~,ord] = sort(mod(atan2(vy-mean(vy),vx-mean(vx)),2*pi));
% vlist_temp = vlist_temp(ord);

C_temp = oldC;

% update the list of the two new cells
V_ref = find(vlist_temp==V_count+1);
ord = [V_ref:z+2,1:V_ref-1];
vlist_new = vlist_temp(ord);
% vlist_new = [vlist_new,vlist_new(1)];
vd = find(vlist_new==V_count+2);
C_temp{cell_id} = vlist_new(1:vd);


V_ref = find(vlist_temp==V_count+2);
ord = [V_ref:z+2,1:V_ref-1];
vlist_new = vlist_temp(ord);
% vlist_new = [vlist_new,vlist_new(1)];
vd = find(vlist_new==V_count+1);
C_temp{Ni+1} = vlist_new(1:vd);

% update the list of the neighboring cells, the edge on the boundary would have only one
% neighboring cell 
if ~isempty(cneib1)
    vid_temp = oldC{cneib1};
    zz = length(vid_temp);

    vid1 = find(vid_temp==v1);
    vid2 = find(vid_temp==v2);
    vid_min = min([vid1,vid2]);
    vid_max = max([vid1,vid2]);
    if vid_min == 1 && vid_max == zz
        vid_temp = [V_count+1,vid_temp];
    elseif vid_min == 1 && vid_max == 2
            vid_temp = [vid_temp(1),V_count+1,vid_temp(2:zz)];
    elseif vid_min ~= 1
        vid_temp = [vid_temp(1:vid_min),V_count+1,vid_temp(vid_min+1:zz)];
    end

    C_temp{cneib1} = vid_temp;
end

if ~isempty(cneib2)
    vid_temp = oldC{cneib2};
    zz = length(vid_temp);

    vid3 = find(vid_temp==v3);
    vid4 = find(vid_temp==v4);
    vid_min = min([vid3,vid4]);
    vid_max = max([vid3,vid4]);
    if vid_min == 1 && vid_max == zz
        vid_temp = [V_count+2,vid_temp];
    elseif vid_min == 1 && vid_max == 2
            vid_temp = [vid_temp(1),V_count+2,vid_temp(2:zz)];
    elseif vid_min ~= 1
        vid_temp = [vid_temp(1:vid_min),V_count+2,vid_temp(vid_min+1:zz)];
    end

    C_temp{cneib2} = vid_temp;
end

inter_points = [intX,intY];
V_temp = [oldV;inter_points];

C_temp = reshape(C_temp,[],1) ;

newV = V_temp;
newC = C_temp;
