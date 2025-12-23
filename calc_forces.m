 function [Cforce,pe,ar] = calc_forces(cx,cy,V,adjCV,newC,K_A,K_P,A0,P0)

[pe, ar] = get_cell_shape(V,newC);
% [pe, ar, ~, ~, ~] = cell_shape_info(V,newC);
Npts = length(newC);
box_size = sqrt(Npts);


% partial derivative of vertex position with respect to cell center position;
dH_mat = zeros(2*2*Npts,2*Npts);
% dH_mat = sparse(dH_mat);
for i = 1:2*Npts
	ind = find(adjCV(:,i)>0);

	% Trio of cell centers define a vertex, which is the circumcenter
	Ax = cx(ind(1));
	Ay = cy(ind(1));
	Bx = cx(ind(2));
	By = cy(ind(2));
	Bx = Ax + pbc_dist(Bx-Ax,box_size);
	By = Ay + pbc_dist(By-Ay,box_size);
	Cx = cx(ind(3));
	Cy = cy(ind(3));
	Cx = Ax + pbc_dist(Cx-Ax,box_size);
	Cy = Ay + pbc_dist(Cy-Ay,box_size);


	% fill dH_matrix
	[HxAx, HxAy, HxBx, HxBy, HxCx, HxCy,...
	HyAx, HyAy, HyBx, HyBy, HyCx, HyCy] = get_dH(Ax,Ay,Bx,By,Cx,Cy);
	dH_mat(i,ind(1)) = HxAx;
	dH_mat(i,ind(1)+Npts) = HxAy;
	dH_mat(i,ind(2)) = HxBx;
	dH_mat(i,ind(2)+Npts) = HxBy;
	dH_mat(i,ind(3)) = HxCx;
	dH_mat(i,ind(3)+Npts) = HxCy;
	dH_mat(i+2*Npts,ind(1)) = HyAx;
	dH_mat(i+2*Npts,ind(1)+Npts) = HyAy;
	dH_mat(i+2*Npts,ind(2)) = HyBx;
	dH_mat(i+2*Npts,ind(2)+Npts) = HyBy;
	dH_mat(i+2*Npts,ind(3)) = HyCx;
	dH_mat(i+2*Npts,ind(3)+Npts) = HyCy;
end

% partial derivative of cell energies with respect to vertex positions;
dE_H_mat = zeros(Npts,2*2*Npts);
% dE_H_mat = sparse(dE_H_mat);
for i = 1:Npts % loop over cells
    neibs = newC{i};
    z = length(neibs);
    % make sure all vertices are within the same PBC image
    Vx = cx(i) + pbc_dist(V(neibs,1)-cx(i),box_size);
    Vy = cy(i) + pbc_dist(V(neibs,2)-cy(i),box_size);
    [~,ord] = sort(mod(atan2(Vy-cy(i),Vx-cx(i)),2*pi));
    neibs = neibs(ord);
    Vx = Vx(ord);
    Vy = Vy(ord);
    
    % make the neigbor list cyclic
    neibs = [neibs, neibs(1), neibs(2)];
    Vx = [Vx; Vx(1); Vx(2)];
    Vy = [Vy; Vy(1); Vy(2)];
    for n = 2:z+1
       Ehx = K_A(i)*(ar(i)-A0(i))*(Vy(n+1)-Vy(n-1)) + ...
             2*K_P(i)*(pe(i)-P0(i))*...
             (Vx(n)-Vx(n-1))/sqrt((Vx(n-1)-Vx(n))^2+(Vy(n-1)-Vy(n))^2) + ...
             2*K_P(i)*(pe(i)-P0(i))*...
             (Vx(n)-Vx(n+1))/sqrt((Vx(n+1)-Vx(n))^2+(Vy(n+1)-Vy(n))^2);
       Ehy = K_A(i)*(ar(i)-A0(i))*(Vx(n-1)-Vx(n+1)) + ...
             2*K_P(i)*(pe(i)-P0(i))*...
             (Vy(n)-Vy(n-1))/sqrt((Vx(n-1)-Vx(n))^2+(Vy(n-1)-Vy(n))^2) + ...
             2*K_P(i)*(pe(i)-P0(i))*...
             (Vy(n)-Vy(n+1))/sqrt((Vx(n+1)-Vx(n))^2+(Vy(n+1)-Vy(n))^2);
       dE_H_mat(i,neibs(n)) = Ehx;
       dE_H_mat(i,neibs(n)+2*Npts) = Ehy;
    end
end

Cforce = zeros(Npts,2);
dH_mat_x = dH_mat(:,1:Npts);
dH_mat_y = dH_mat(:,Npts+1:2*Npts);
Cforce(:,1) = -sum(dE_H_mat*dH_mat_x,1);
Cforce(:,2) = -sum(dE_H_mat*dH_mat_y,1);
end