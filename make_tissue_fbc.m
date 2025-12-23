function  [newV, newC] = make_tissue_fbc(centers, Ni, V, cell_list)

Npts = size(cell_list,1);
box_size = sqrt(Npts);
cx0 = box_size/2;
cy0 = box_size/2;
cx = centers(:,1);
cy = centers(:,2);
L_list = sqrt((cx-cx0).^2 + (cy-cy0).^2);

% choose the number of inner cells
[~,idx]=sort(L_list);
in_cells = idx(1:Ni);
% out_cells = setdiff(1:Npts,in_cells);
% No = Npts - ii;

clist = cell(Ni,1);
vertex_position = [];
for i = 1:Ni
    in_cid = in_cells(i);
    x0 = cx(in_cid);
    y0 = cy(in_cid);
    vlist = cell_list{in_cid};
    vx_temp = V(vlist,1);
    vy_temp = V(vlist,2);
    vx = x0 + pbc_dist(vx_temp-x0,box_size);
    vy = y0 + pbc_dist(vy_temp-y0,box_size);
    V_cell = [vx,vy];
    Z = length(vlist);
    clist_temp = zeros(1,Z);
    for j = 1:Z
        v_refer = V_cell(j,:);
        if i == 1
            vertex_position = [vertex_position; v_refer];
            clist_temp(j) = j;
        elseif i > 1
            kk = length(vertex_position);
            vertex_position_new = zeros(kk,2);
            for k = 1:kk
                v2_temp = vertex_position(k,:);
                vertex_position_new(k,:) = v2_temp;
            end
            
            diff_temp = v_refer - vertex_position_new;
            distance_list = sqrt(sum(diff_temp.^2,2));
            
            if min(distance_list)<1e-10 % no overlap
                [~, min_pos] = min(distance_list);
                clist_temp(j) = min_pos;
            else
                vertex_position = [vertex_position; v_refer];
                clist_temp(j) = size(vertex_position, 1);
            end
        end
        clist{i} = clist_temp;
    end
end

newV = vertex_position;
newC = clist;

end