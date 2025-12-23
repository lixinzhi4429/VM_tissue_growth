function [cx,cy,P0,i_relax,max_force,final_energy] = spv_gs_FIRE(mean_p0,std_p0,seed_number,Npts,force_tol,n_relax,cx,cy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set Random Seed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = RandStream('mcg16807','Seed',seed_number);
RandStream.setGlobalStream(s);


box_size = sqrt(Npts); % Linear size of the bounding box, should be sqrt(Npts)

if cx == 0 && cy == 0
    cx = rand(Npts,1);
    cy = rand(Npts,1);
    cx = box_size*cx;
    cy = box_size*cy;
end
cx = cx(:);
cy = cy(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_A = ones(Npts,1); % set modulus for perimeter term to be one
K_P = ones(Npts,1); % set modulus for area term to be one
A0 = ones(Npts,1); % A0 doesn't MATTER!!!
P0 = std_p0*randn(Npts,1) + mean_p0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FIRE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fire parameters
Nmin = 5;
finc = 1.1;
fdec = 0.5;
f_alpha = 0.99;
alpha_start = 0.1;
dt_max = 0.1;
dt_init = 0.01;
fire_mass = 4;
% force_tol = 1e-10;
% n_relax = 2e4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = dt_init;
vx = zeros(Npts,1);
vy = zeros(Npts,1);
alpha = alpha_start;
max_force = 1e8;
i_relax = 0;
Nneg = 0;

while max_force > force_tol && i_relax < n_relax 

    [V, newC, ~, adjCV, ~] = make_voronoi(cx/box_size, cy/box_size);
    [Cforce,~,~]  = calc_forces(cx,cy,V,adjCV,newC,K_A,K_P,A0,P0); 
    fx = Cforce(:,1);
    fy = Cforce(:,2);
    
    fmag = sqrt(fx.^2+fy.^2);
    fxhat = fx./fmag;
    fyhat = fy./fmag;
    fxhat(isnan(fxhat)) = 0;
    fyhat(isnan(fyhat)) = 0;
    Pfire = sum(vx.*fx+vy.*fy);
    




    if Pfire > 0 && Nneg > Nmin
        dt = min(dt*finc,dt_max);
        alpha = alpha*f_alpha;
        vmag = sqrt(vx.^2+vy.^2);
        vx = (1-alpha)*vx + alpha* fxhat .* vmag;
        vy = (1-alpha)*vy + alpha* fyhat .* vmag;
    end
    if Pfire <= 0
        dt = dt*fdec;
        alpha = alpha_start;
        vx = zeros(Npts,1);
        vy = zeros(Npts,1);
        Nneg = 0;
    end
    Nneg = Nneg+1;
    vx = vx + fx*dt/fire_mass;
    vy = vy + fy*dt/fire_mass;
    cx = mod(cx + vx*dt,box_size);
    cy = mod(cy + vy*dt,box_size);
    

    max_force = mean(sqrt(sum(Cforce.^2,2)));
    if mod(i_relax,100) == 0
%         e0 = sum(K_P.*(pe-P0).^2+K_A.*(ar-A0).^2)/Npts;
        disp(['step =',num2str(i_relax)]);
        disp(['current max force =',num2str(max_force,'%.5e\n')]);
%         disp(['current energy =',num2str(e0,'%.5e\n')]);

    end
    i_relax = i_relax + 1;

end


% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % LBFGS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % LBFGS parameters
% invH = eye(2*Npts);
% dt0 = 1;
% i_relax = 0;
% max_force = 1e8;
% % force_tol = 1e-10;
% % n_relax = 2e4;
% 
% 
% 
% rho_list = sort([(1/2).^(0:2:15)]);
% n_lnsrch = length(rho_list);
% de = zeros(n_lnsrch,2);
% de(:,1) = dt0*rho_list;
% while max_force > force_tol && i_relax < n_relax 
%     [V, newC, ~, adjCV, ~] = make_voronoi(cx/box_size, cy/box_size);
%     [Cforce,pe,ar]  = calc_forces(cx,cy,V,adjCV,newC,K_A,K_P,A0,P0);
%     ff = [Cforce(:,1);Cforce(:,2)];
%     e0 = sum(K_P.*(pe-P0).^2+K_A.*(ar-A0).^2)/Npts;
%  
%     for j = 1:n_lnsrch
%         
%         dc = de(j,1)*invH*ff;
%         cxi = mod(cx + dc(1:Npts),box_size);
%         cyi = mod(cy + dc((Npts+1):2*Npts),box_size);
%         [V, newC, ~, ~, ~] = make_voronoi(cxi/box_size, cyi/box_size);
%         [pe, ar] = get_cell_shape(V,newC);
%         enew = sum(K_P.*(pe-P0).^2+K_A.*(ar-A0).^2)/Npts;
%         de(j,2) = enew-e0;
%     end
%     
%     [~,b]=min(de(:,2));
%     dt = de(b,1);
%     
%     dc = dt*invH*ff;
%     cx = mod(cx + dc(1:Npts),box_size);
%     cy = mod(cy + dc(Npts+1:2*Npts),box_size);
% 
%     oldff = ff;
%     [V, newC, ~, adjCV, ~] = make_voronoi(cx/box_size, cy/box_size);
%     [Cforce,~,~]  = calc_forces(cx,cy,V,adjCV,newC,K_A,K_P,A0,P0);
%     ff = [Cforce(:,1);Cforce(:,2)];
%     dff = ff-oldff;
%     
%     AA = -(dc*dc')/(dc'*dff);
%     kk = (dff'*(invH*dff));
%     BB = - (invH*dff)*(invH*dff)'/kk;
%     uu = -dc/(dc'*dff) + (invH*dff)/kk;
%     CC = kk*(uu*uu');
%     invH = invH + AA+BB+CC;
%     disp(invH(1,1));
%     max_force = mean(sqrt(sum(Cforce.^2,2)));
%     if mod(i_relax,20) == 0
%         disp(['step =',num2str(i_relax)]);
%         disp(['current max force =',num2str(max_force,'%.5e\n')]);
%         disp(['current energy =',num2str(e0,'%.5e\n')]);
% 
%     end
% %     disp(sum((K_A.*(ar-A0).^2))/Npts);
%     i_relax = i_relax + 1;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Save Linear Response Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% final_track_file = [output_path,'/XY_npts_',num2str(Npts),'_seed_',num2str(seed_number),'_p0_',num2str(p0_mean),'_',num2str(p0_width),'.txt'];
% final_track_file = [output_path,'/',output_filename];
% final_track_file = output_filename;

% dlmwrite(final_track_file,[cx, cy, K_A, K_P, A0, P0],'precision', 16);
[~,pe,ar]  = calc_forces(cx,cy,V,adjCV,newC,K_A,K_P,A0,P0); 
final_energy = sum(K_P.*(pe-P0).^2+K_A.*(ar-A0).^2);

cx = cx/box_size;
cy = cy/box_size;

end