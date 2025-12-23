clear all;
close all;
% numcores = feature('numcores');
% parpool(numcores);
% parpool(40);

% filename too long will cause error, like section %% does not work

Npts = 400;
box_size = sqrt(Npts); % square box size of the system
L = box_size;

mean_p0 = 3.75;
std_p0 = 0;
seed_number = 1;
force_tol = 1e-10;

[cx,cy,p0_list,~,~,~] = spv_gs_FIRE(mean_p0,std_p0,seed_number,Npts,force_tol,500,0,0);

cx = cx*box_size;
cy = cy*box_size;

[V, cell_list, ~, ~, ~] = make_voronoi_ccw(cx/box_size, cy/box_size);

cm = jet(Npts);
figure(1);
clf
for i=1:Npts
    X = V(cell_list{i},:);
    Z = size(X,1);
    Ximg = zeros(Z,2);
    Ximg(1,:) = X(1,:);

    Ximg(2:Z,1) = X(1,1) + pbc_dist(X(2:Z,1)-X(1,1),box_size);
    Ximg(2:Z,2) = X(1,2) + pbc_dist(X(2:Z,2)-X(1,2),box_size);

    LL = plot([Ximg(:,1);Ximg(1,1)],[Ximg(:,2);Ximg(1,2)],'k-','linewidth',2);
    hold on;
    h = patch([Ximg(:,1);Ximg(1,1)],[Ximg(:,2);Ximg(1,2)],'w','facecolor',cm(i,:));
    hold on;

end
%% make initial state under fbc
centers = [cx, cy];
N0 = 1;
[newV, newC] = make_tissue_fbc(centers, N0, V, cell_list);
V_count = size(newV,1);
% plot the cut structure
cm = jet(N0);
cell_centers = zeros(N0,2);
figure(2);
clf;
for i = 1:N0
    z = length(newC{i});
    vx = newV(newC{i},1);
    vy = newV(newC{i},2);
    cell_centers(i,:) = mean([vx,vy],1);
    plot(vx([1:end,1]),vy([1:end,1]),'k-','linewidth',2);
    hold on;
    patch(vx([1:end,1]),vy([1:end,1]),'w','facecolor','b');
    hold on;
    txt = num2str(i);
    text(cell_centers(i,1)+0.05,cell_centers(i,2),txt,'fontsize',20);
    hold on;

end

plot(cell_centers(:,1),cell_centers(:,2),'k.', 'markersize',20,'linewidth',2);

for j = 1:V_count
    txt = num2str(j);
    text(newV(j,1)+0.05,newV(j,2),txt,'fontsize',20,'color','green');
    hold on;
end

plot(newV(:,1),newV(:,2),'ro', 'markersize',5,'linewidth',2);
hold on;

pbaspect([1,1,1]);
height=1; % width/golden ratio
width=1;
% scale=800; % 600 is 3.13 inches
% xpos=50;
% ypos=500;
% box off;
% xlim([-1,6]);
% ylim([-1,6]);
% set(gca, 'YTick', []);
% set(gca, 'XTick', []);
% % set(gca,'Visible','off');
% % set(gcf,'Position',[xpos ypos scale*width scale*height]);
% set(gca,'FontName','Times New Roman');
% set(gca,'FontSize',50);
% set(gca,'LineWidth',2);
% set(gca, 'Color', 'none'); % Sets axes background
% set(0,'DefaultFigureColor',[1 1 1])
set(gca, 'Color', 'none'); % Sets axes background
set(gcf,'color','w');

% fn = ['snapshot_N_',num2str(Npts),'_seed_',num2str(seed),'_p0_',num2str(mean_p0),'_std_',num2str(std_p0),'.pdf'];
ax = gca;
% exportgraphics(ax,fn,'ContentType','vector','BackgroundColor','none');

[nE,E_list, adjCV, adjVV, adjCC, adjCE, adjVE] = make_adj_fbc(newV, newC);

figure(2);
% cids = find(sum(adjCC,1)<=3);
% cc = length(cids);
% for c2 = 1:cc
%     cid = cids(c2);
%     z = length(newC{cid});
%     vx = newV(newC{cid},1);
%     vy = newV(newC{cid},2);
%     plot(vx([1:end,1]),vy([1:end,1]),'k-','linewidth',2);
%     hold on;
%     patch(vx([1:end,1]),vy([1:end,1]),'w','facecolor','g');
% end

eids = find(sum(adjCE,1)==1);
ee = length(eids);
for j = 1:ee
    eid = eids(j);
    Ximg = newV(E_list(eid,:),:);
    plot([Ximg(1,1);Ximg(2,1)],[Ximg(1,2);Ximg(2,2)],'r-','MarkerSize',2,'linewidth', 6);
    hold on;
end

vids = find(sum(adjCV,1)<=2);
vv = length(vids);
for j = 1:vv
    vid = vids(j);
    vx = newV(vid,1);
    vy = newV(vid,2);
    plot(vx,vy,'ro','MarkerSize',10,'linewidth',4);
    hold on;
end

boundary_vertices = vids;

pbaspect([1,1,1]);
height=1; % width/golden ratio
width=1;
% scale=800; % 600 is 3.13 inches
% xpos=5
% ypos=500;
box on;
% xlim([0,10]);ai
% ylim([0,10]);
% set(gca,'XTick',0:2:10);
% set(gca,'YTick',0:2:10);
% xticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'});
% yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'});
% % set(gca,'Visible','off');
% % set(gcf,'Position',[xpos ypos scale*width scale*height]);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',50);
set(gca,'LineWidth',2);
% set(gca, 'Color', 'none'); % Sets axes background
% set(0,'DefaultFigureColor',[1 1 1]);
set(gca, 'Color', 'none'); % Sets axes background
set(gcf,'color','w');
% set(gcf, 'Position', get(0, 'Screensize'));% make figure full screen

fn = ['snapshot_N_',num2str(N0),'_Ni_',num2str(N0),'_seed_',num2str(seed_number),'_p0_',num2str(mean_p0),'_std_',num2str(std_p0),'_pickup_inncer_cells_boundary.pdf'];
ax = figure(2);
% exportgraphics(ax,fn,'ContentType','vector','BackgroundColor','none');

[perim_list, area_list] = get_cell_shape_fbc(newV,newC);
%% relax the state

V_list_int = newV;
Cell_list_int = newC;

% the large cell include all boundary vertices and acts as medium,
% try the large cell separately in the simulations
% make sure cell_large in clockwise direction,
% be careful, cell large is in clockwise direction, different from other
% cells, the tricky thing.

cell_large_int = boundary_vertices;
Z = length(cell_large_int);
for i = 1:Z
    Vx = V_list_int(cell_large_int,1);
    Vy = V_list_int(cell_large_int,2);
    [~,ord] = sort(mod(atan2(Vy-mean(Vy),Vx-mean(Vx)),2*pi));
    cell_large_int = cell_large_int(ord);
end

cell_large_int = flip(cell_large_int);


KP = 0.01;
KA_vals = [1:2:10,10];
KA_list = repmat(KA_vals,1,6);
kk = length(KA_list);
k_vals = [10,20,50,100,150,200];
k_list = repmat(k_vals,6,1);
k_list = k_list(:)';

dt_list = [0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.05,0.01,0.01,0.01,0.01,0.01,0.05,0.01,0.01,0.01,0.01,0.01]/2;
G2_check = 5;

sigma = 0;
tau = 0.1;
%%
Dr = 1;
v0 = 0.1;
V_count_temp = size(V_list_int,1);
theta_int = 2*pi*rand(N0,1);

KA = 1;
k = 100;
dt = 0.01;
KA_large = 0; % make sure the large medium exert no pressure on inner cells
KP_large = 0.01;
A0_initital = 1;
A0 = A0_initital;
A0_large = 1;

N0 = size(Cell_list_int,1);
K_P = KP*ones(N0+1,1);
K_A = KA*ones(N0+1,1);
K_P(N0+1) = KP_large;
K_A(N0+1) = KA_large;
% A0_list = zeros(Ni+1,1);
% A0_list(1:end-1) = A0;
% A0_list(end) = A0_large; %Npts - sum(A0_list);

% A0_list = [0.75:0.05:1.2]';

%     A0_list_int = linspace(2.4,3.2,N0)';
A0_list_int = 1;
A0_list_int(N0+1,1) = A0_large;

Aids = 1:N0;
Bids = N0 + 1;

lambda_0 = -0.02;
lambda_AA = lambda_0;
lambda_AB = lambda_0;

[nE,E_list, ~, ~,~,~,adjCE,~] = make_adj_fbc_large_cell(V_list_int, Cell_list_int,cell_large_int);
lambda_list = zeros(nE,1);
for eid2 = 1:nE
    cneib = find(adjCE(:,eid2));
    c1 = cneib(1);
    c22 = cneib(2);
    if sum(c1==Aids)>0 && sum(c22==Aids)>0
        lambda_list(eid2) = lambda_AA;
    elseif sum(c1==Aids)>0 && sum(c22==Bids)>0 || sum(c1==Bids)>0 && sum(c22==Aids)>0
        lambda_list(eid2) = lambda_AA;
    end

end

force_tol = 1e-12;
i_relax = 0;
max_force = 1e8;

V_list = V_list_int;
Cell_list = Cell_list_int;
cell_large = cell_large_int;
A0_list = A0_list_int;
theta = theta_int;

dt_relax = 0.02;
while i_relax <= 20000 && max_force > force_tol

    [Vertex_Force, ~, Cell_area] = get_VM_force_fbc_tension_large_cell(Npts,V_list, Cell_list,cell_large,K_P,K_A,A0_list,lambda_list,E_list);
    Vx = V_list(:,1);
    Vy = V_list(:,2);
    fx = Vertex_Force(:,1);
    fy = Vertex_Force(:,2);

    Vx = Vx + fx*dt_relax;
    Vy = Vy + fy*dt_relax;

    V_list = [Vx,Vy];

    max_force = max(sqrt(sum(Vertex_Force.^2,2)));
    if mod(i_relax,100) == 0
        disp(['step =',num2str(i_relax)]);
        disp(['current max force =',num2str(max_force,'%.5e\n')]);
    end

    i_relax = i_relax + 1;

end

[~, area] = get_cell_shape_fbc(V_list,Cell_list);
% try to include varying A0 and cell division
% update newC to include the large medium cell.

% change A0 growth rule to incorporate contact inhibition of proliferation

output_dir = './test_data_cell_polar/';
n_relax = 100000;
force_tol = 1e-12;

div_rate = 1;
gamma = 1;

std_g = 0;
mean_g = 0.5;

max_force = 1e8;
T1_step = 1;
T1_threshold = 0.01;

% when cell area reaches V_sizer, enter G2 phase
G1_sizer = 2;

% [~,~,~,adjVV,~,~,~] = make_adj_fbc(V_list, Cell_list);
% adj_lambda = lambda_0*adjVV;
% % logical(adj_lambda), adjVV should be equal
% isequal(logical(adj_lambda), adjVV);


% cells divide after tau steps in G2 phase
% specify tau using relax step, real timer t = tau*dt
G2_timer = G2_check/dt;

division_fn = ['./area_data_cell_polar/division_Ni_',num2str(N0),'_seed_',num2str(seed_number),'_th_',num2str(T1_threshold),'_sizer_',num2str(G1_sizer),'_timer_',num2str(G2_timer),'_gamma_',num2str(gamma),'_nsteps_',num2str(n_relax),'_A0_1','_dt_',num2str(dt),'_mean_G_',num2str(mean_g),'_k_',num2str(k),'_new_A0_pressure_KA_large_',num2str(KA_large),'_KP_',num2str(KP),'_KA_',num2str(KA),'_lambda0_',num2str(lambda_0),'_sigma_',num2str(sigma),'_tau_',num2str(tau),'_Dr_',num2str(Dr),'_v0_',num2str(v0),'_div_large_round_bc.txt'];
T1_fn = ['./area_data_cell_polar/T1_Ni_',num2str(N0),'_seed_',num2str(seed_number),'_th_',num2str(T1_threshold),'_sizer_',num2str(G1_sizer),'_timer_',num2str(G2_timer),'_gamma_',num2str(gamma),'_nsteps_',num2str(n_relax),'_A0_1','_dt_',num2str(dt),'_mean_G_',num2str(mean_g),'_k_',num2str(k),'_new_A0_pressure_KA_large_',num2str(KA_large),'_KP_',num2str(KP),'_KA_',num2str(KA),'_lambda0_',num2str(lambda_0),'_sigma_',num2str(sigma),'_tau_',num2str(tau),'_Dr_',num2str(Dr),'_v0_',num2str(v0),'_div_large_round_bc.txt'];
area_fn = ['./area_data_cell_polar/area_Ni_',num2str(N0),'_seed_',num2str(seed_number),'_th_',num2str(T1_threshold),'_sizer_',num2str(G1_sizer),'_timer_',num2str(G2_timer),'_gamma_',num2str(gamma),'_nsteps_',num2str(n_relax),'_A0_1','_dt_',num2str(dt),'_mean_G_',num2str(mean_g),'_k_',num2str(k),'_new_A0_pressure_KA_large_',num2str(KA_large),'_KP_',num2str(KP),'_KA_',num2str(KA),'_lambda0_',num2str(lambda_0),'_sigma_',num2str(sigma),'_tau_',num2str(tau),'_Dr_',num2str(Dr),'_v0_',num2str(v0),'_div_large_round_bc.txt'];


G2_clock = zeros(N0,1);
cells_G2 = find(area>G1_sizer);
ll = length(cells_G2);

% set RandStream to make sure each step use same random numbers.
s = RandStream('mcg16807','Seed',seed_number);
RandStream.setGlobalStream(s);
G2_clock_initial = randperm(G2_timer,ll);
G2_clock(cells_G2) = G2_clock_initial';

g_list = mean_g*ones(N0,1);

i_relax = 0;
div_count = 0;
T1_count = 0;

[~,~,~,adjVV,~,~,~] = make_adj_fbc(V_list, Cell_list);
adj_lambda = lambda_0*adjVV;

% try
while i_relax <= n_relax

    [~, ar] = get_cell_shape_fbc(V_list,Cell_list);

    % pick up cells entering G2 phase，start G2 clock
    G2_cells = find(ar>G1_sizer | G2_clock>0);
    G2_clock(G2_cells) = G2_clock(G2_cells) + 1;

    V_count_old = size(V_list,1);
    area_temp = ar;
    divide_cells = find(G2_clock>=G2_timer);

    cell_large_temp = flip(cell_large);
    vx_large = V_list(cell_large,1);
    vy_large = V_list(cell_large,2);
    [geom_t,~, cpmo_t] = polygeom(vx_large, vy_large);

    tissue_x_cen = geom_t(2);
    tissue_y_cen = geom_t(3);
    I1_t = cpmo_t(1);
    angle1_t = cpmo_t(2);
    I2_t = cpmo_t(3);
    angle2_t = cpmo_t(4);

    tissue_center = mean(V_list,1);

    Cell_list_temp = Cell_list;
    V_list_temp = V_list;

    if mod(i_relax,div_rate) == 0 && ~isempty(divide_cells) %&& i_relax<6000
        for i = 1:length(divide_cells)
            N_old = size(Cell_list,1);
            V_count_old = size(V_list,1);
            cid = divide_cells(i);

            div_center = mean(V_list_temp(Cell_list_temp{cid},:),1);
            pressure = -KA*(area_temp(cid) - A0_list(cid));

            vx_temp = V_list_temp(Cell_list_temp{cid},1);
            vy_temp = V_list_temp(Cell_list_temp{cid},2);
            [geom,~, cpmo] = polygeom(vx_temp, vy_temp);
            x_cen = geom(2);
            y_cen = geom(3);
            I1 = cpmo(1);
            angle1 = cpmo(2);
            I2 = cpmo(3);
            angle2 = cpmo(4);

            [V_list, Cell_list,vertices_split] = make_cell_division_new_lambda(cid, V_list, Cell_list);
            div_count = div_count + 1;

            G2_clock(cid) = 0;
            G2_clock(N_old + 1,1) = 0;

            % divide A0 based on the real area of new cells

            [~, area_temp] = get_cell_shape_fbc(V_list,Cell_list);

            area_1 = area_temp(cid);
            area_2 = area_temp(N_old+1);
            area_total = (area_1+area_2);

            %                     writematrix([i_relax,cid,N_old+1,div_count,area_total,area_1,area_2],division_fn,'WriteMode','append');
            % writematrix([i_relax,cid,N_old+1,div_count,area_1,area_2,div_center,I1,angle1,I2,angle2,tissue_center,I1_t,angle1_t,I2_t,angle2_t],division_fn,'WriteMode','append');



            A0_list(cid) = area_1 + pressure/KA;
            A0_list(N_old+1) = area_2 + pressure/KA;

            g_temp = g_list(cid);
            g_list(cid,1) = g_temp;
            g_list(N_old+1,1) = g_temp;

            % theta(V_count_old+1) = 2*pi;
            % theta(V_count_old+2) = 2*pi;

            theta_rand = 2*pi*rand(N_old+1,1);
            theta(cid) = theta_rand(cid);
            theta(N_old+1,:) = theta_rand(N_old+1);

            % P_vec(cid,:) = [cos(theta_rand(cid)),sin(theta_rand(cid))];
            % P_vec(N_old+1,:) = [cos(theta_rand(N_old+1)),sin(theta_rand(N_old+1))];

            % update tension adjacency matrix
            [~,~,adjCV,adjVV,~,~,~] = make_adj_fbc(V_list, Cell_list);
            adj_lambda_old = adj_lambda;

            % disconnect old edges between which new vertices are added
            v1 = vertices_split(1);
            v2 = vertices_split(2);
            v3 = vertices_split(3);
            v4 = vertices_split(4);

            adj_lambda_old(v1,v2) = 0;
            adj_lambda_old(v2,v1) = 0;
            adj_lambda_old(v3,v4) = 0;
            adj_lambda_old(v4,v3) = 0;

            V_count = size(V_list,1);
            adj_lambda_temp = lambda_0*adjVV;
            adj_lambda_temp(1:V_count_old,1:V_count_old) = adj_lambda_old.*adjVV(1:V_count_old,1:V_count_old);
            adj_lambda = adj_lambda_temp;

            %             if ~isequal(logical(adj_lambda), adjVV)
            %                 disp([i_relax, cid]);
            %                 break;
            %             end
            % update cell_large list
            % free boundary adjacency matrices, boundary vertices have two
            % neighboring cells
            v_bound = sum(adjCV(:,V_count-1:V_count));

            if v_bound(1) == 2 && v_bound(2) ~= 2
                vid_bound = V_count -1;
            elseif v_bound(1) ~= 2 && v_bound(2) == 2
                vid_bound = V_count;
            elseif v_bound(1) == 2 && v_bound(2) == 2
                vid_bound = [V_count-1,V_count];
            elseif sum(v_bound==2) == 0
                % new vertices not on the boundary
                % no need to update cell_large
                continue;
            end

            for vv2 = 1:length(vid_bound)
                vid = vid_bound(vv2);
                vids_temp = find(adjVV(:,vid));
                vids = intersect(vids_temp,cell_large);

                id1 = find(cell_large==vids(1));
                id2 = find(cell_large==vids(2));
                id = min([id1,id2]);
                id_temp = max([id1,id2]);
                if id == 1 && id_temp == 2
                    cell_large = [cell_large(1),vid,cell_large(2:end)];
                elseif id == 1 && id_temp == length(cell_large)
                    cell_large = [vid, cell_large];
                elseif id ~= 1
                    cell_large = [cell_large(1:id), vid,cell_large(id+1:end)];
                end
            end
        end

    end


    [~, ar] = get_cell_shape_fbc(V_list,Cell_list);

    Ni = size(Cell_list,1);
    K_P = KP*ones(Ni+1,1);
    K_A = KA*ones(Ni+1,1);
    K_P(Ni+1) = KP_large;
    K_A(Ni+1) = KA_large;

    A0_list(1:Ni) = A0_list(1:Ni) + g_list.*exp(-k.*(ar - A0_list(1:Ni)).^2)*dt;
    A0_list(Ni+1) = 1; % large cell target area, could be any value because KA_large = 0

    Aids = 1:Ni;
    Bids = Ni + 1;

    [nE,E_list,lambda_list,adjCV,~,~,~] = make_adj_fbc_lambda(V_list, Cell_list,adj_lambda);


    [Vertex_Force, ~, Cell_area] = get_VM_force_fbc_tension_large_cell(Npts,V_list, Cell_list,cell_large,K_P,K_A,A0_list,lambda_list,E_list);
    Vx = V_list(:,1);
    Vy = V_list(:,2);
    fx = Vertex_Force(:,1);
    fy = Vertex_Force(:,2);

    % Vx = Vx + gamma*fx*dt;
    % Vy = Vy + gamma*fy*dt;

    % V_count = size(V_list,1);
    % Dr_list = Dr*ones(V_count,1);
    % v0_list = v0*ones(V_count,1);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % Polarization Dynamics (VECTORIZED)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % theta = theta + randn(V_count,1).*sqrt(2*dt.*Dr_list);

    % fvx = gamma*fx + v0_list.*cos(theta);
    % fvy = gamma*fy + v0_list.*sin(theta);


    Dr_list = Dr*ones(Ni,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Polarization Dynamics (VECTORIZED)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = theta + randn(Ni,1).*sqrt(2*dt.*Dr_list);

    %%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Cell Dynamics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_vec = [cos(theta), sin(theta)];
    % get vertex polarity
    V_count = size(V_list,1);
    PV_vec = zeros(V_count,2);
    for j = 1:V_count
        cids = adjCV(:,j);
        PV_vec(j,:) = mean(P_vec(cids,:),1);
    end

    fvx = gamma*fx + v0*PV_vec(:,1);
    fvy = gamma*fy + v0*PV_vec(:,2);

    Vx = Vx + gamma*fvx*dt;
    Vy = Vy + gamma*fvy*dt;

    %     % implement round fixed boundary condition
    %     Vx0 = box_size/2;
    %     Vy0 = box_size/2;
    %
    %     dx = Vx - Vx0;
    %     dy = Vy - Vy0;
    %     dr = sqrt(dx.^2+dy.^2);
    %     R = box_size/2;
    %     vids = find(dr>R);
    %     vv = length(vids);
    %     for v2 = 1:vv
    %         vid = vids(v2);
    %         dx_temp = dx(vid);
    %         dy_temp = dy(vid);
    %
    %         if dx_temp ~= 0
    %             theta = mod(atan2(dy_temp,dx_temp),2*pi);
    %             Vx(vid) = R*cos(theta) + R;
    %             Vy(vid) = R*sin(theta) + R;
    %         elseif dx_temp == 0 && dy_temp > R
    %             Vx(vid) = R;
    %             Vy(vid) = 2*R;
    %         elseif dx_temp == 0 && dy_temp < 0
    %             Vx(vid) = R;
    %             Vy(vid) = 0;
    %         end
    %
    %     end

    V_list = [Vx,Vy];

    [~,E_list,edge_length_list,~,adjVV2,~,~,~,~] = make_adj_fbc_T1(V_list, Cell_list);
    eids = find(edge_length_list<T1_threshold);
    ee = length(eids);
    if ~isempty(eids)
        for e2 = 1:ee
            eid = eids(e2);
            [V_list, Cell_list, cell_large,adj_lambda,T1_count] = make_T1_fbc_large_cell_lambda(V_list, Cell_list, cell_large, eid,adj_lambda,lambda_0,T1_count);
            % writematrix([i_relax,eid,E_list(eid,:),T1_count],T1_fn,'WriteMode','append');
        end
    end
    [pe, ar] = get_cell_shape_fbc(V_list,Cell_list);
    [~,~,~,~,adjVV,~,~] = make_adj_fbc_lambda(V_list, Cell_list,adj_lambda);
    vv = size(adj_lambda,1);
    % adj_lambda = adj_lambda - dt/tau.*(adj_lambda-lambda_0.*adjVV)+ randn(vv).*adjVV.*sqrt(2*dt.*sigma^2/tau);


    [nE,E_list,lambda_list,~,adjVV,~,adjCE] = make_adj_fbc_lambda(V_list, Cell_list,adj_lambda);
    test_equal = isequal(logical(adj_lambda), adjVV);

    lambda_list = lambda_0*ones(nE,1);
    tension_all = (KP*pe'*adjCE + lambda_list')';
    pressure_cell = -KA*(ar - A0_list(1:end-1));

    tension_cell = zeros(Ni,1);
    for c2 = 1:Ni
        eids_temp = find(adjCE(c2,:));
        tension_cell(c2) = mean(tension_all(eids_temp));
    end

    % if Ni <= 20
    %     writematrix([i_relax,Ni,ar'],area_fn,'WriteMode','append');
    % else
    %     writematrix([i_relax,Ni,ar(1:20)'],area_fn,'WriteMode','append');
    % end

    % save configuration information

    vm_info = struct(...
        'Ni',[],...
        'step',[],...
        'V_list',[],...
        'Cell_list',[],...
        'cell_large',[],...
        'cell_shapes',[],...
        'cell_polarity',[],...
        'lambda_matrix',[],...
        'lambda_list',[],...
        'cell_forces',[],...
        'cell_params',[],...
        'div_params',[],...
        'tension_fluc_params',[],...
        'T1_count',[],...
        'G2_clock',[]...
        );


    vm_info.Ni = Ni;
    vm_info.step = i_relax;
    vm_info.V_list = V_list;
    vm_info.Cell_list = Cell_list;
    vm_info.cell_large = cell_large;
    vm_info.cell_shapes = [pe,ar,A0_list(1:Ni),theta];
    vm_info.cell_polarity = theta;
    vm_info.lambda_matrix = adj_lambda;
    vm_info.lambda_list = lambda_list;
    vm_info.cell_forces = [tension_cell,pressure_cell];
    vm_info.cell_params = [KP,KA];
    vm_info.div_params = [dt,div_rate,T1_step,T1_threshold,k,gamma,mean_g,G1_sizer,G2_timer,Dr,v0];
    vm_info.tension_fluc_params = [lambda_0,sigma,tau];
    vm_info.T1_count = T1_count;
    vm_info.G2_clock = G2_clock;


    % if mod(i_relax,10) == 0
    %     all_info_mat = [output_dir,'Ni_',num2str(N0),'_seed_',num2str(seed_number),'_th_',num2str(T1_threshold),'_sizer_',num2str(G1_sizer),'_timer_',num2str(G2_timer),'_gamma_',num2str(gamma),'_A0_1','_dt_',num2str(dt),'_mean_G_',num2str(mean_g),'_k_',num2str(k),'_new_A0_pressure_KA_large_',num2str(KA_large),'_step_',num2str(i_relax),'_KP_',num2str(KP),'_KA_',num2str(KA),'_lambda0_',num2str(lambda_0),'_sigma_',num2str(sigma),'_tau_',num2str(tau),'_Dr_',num2str(Dr),'_v0_',num2str(v0),'_div_large_round_bc.mat'];
    %     parsave(all_info_mat,vm_info);
    %     disp([div_count,test_equal,i_relax]);
    % end

    i_relax = i_relax + 1;
    if mod(i_relax,100) == 0
        %     %         fprintf('%d %d %0.5f %0.5f\n',div_count,i_relax,max(area_list),mean(area_list))
        %     disp([div_count,test_equal,i_relax]);


        growth_rate = g_list.*exp(-k.*(ar - A0_list(1:Ni)).^2);
        [~,ids] = sort(growth_rate);

        c1 = 0;
        c2 = 1;
        mymap = [c1*ones(Ni,1),linspace(0,1,Ni)',c2*ones(Ni,1) ];
        % mymap = jet(N);
        cell_centers = zeros(Ni,2);
        figure(1);
        clf;
        for i = 1:Ni
            z = length(Cell_list{i});
            vx = V_list(Cell_list{i},1);
            vy = V_list(Cell_list{i},2);
            cell_centers(i,:) = mean([vx,vy],1);
            plot(vx([1:end,1]),vy([1:end,1]),'k-','linewidth',2);
            hold on;
            color_id = find(ids==i);
            patch(vx([1:end,1]),vy([1:end,1]),'w','facecolor',mymap(color_id,:));
            hold on;
        end
        ff = 0.1;
        % quiver(cell_centers(:,1),cell_centers(:,2),ff*P_vec(:,1),ff*P_vec(:,2),0,'linewidth',3,'color','black');
        if Ni > 2

            colormap(mymap);
            c = colorbar;
            c.FontSize = 20;

            c.Title.String = "Growth";
            c.FontSize = 30;
            clim([min(growth_rate) max(growth_rate)]);
            pbaspect([1,1,1]);
            % xlim([-25,50]);
            % ylim([-25,50]);
            xlim([-10,25]);
            ylim([-10,25]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'XColor', 'none','YColor','none');
            set(gca, 'Color', 'none'); % Sets axes background
            height=1.0/1.0; % width/golden ratio
            width=1;
            scale=600; % 600 is 3.13 inches
            xpos=50;
            ypos=50;
            set(gcf,'Position',[xpos ypos scale*width scale*height]);
            set(gcf,'color','w');
            box off;
            axis equal
            % %
            %         drawnow
            %         currFrame = getframe(gcf);
            %         writeVideo(video_object,currFrame);

        end

    end
    % close(video_object);
end
