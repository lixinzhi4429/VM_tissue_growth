function [nE,E_list,adjCV, adjVV, adjCC, adjCE, adjVE] = make_adj_fbc(newV, newC)

V_count = size(newV,1);
Ni = size(newC,1);
% ordered Cell-Vertex adjacency matrix 
adjCV = zeros(Ni,V_count);
for i = 1:Ni
    zz = length(newC{i});
    adjCV(i,newC{i}) = 1:zz;
end
adjCV = adjCV > 0;

adjVV = zeros(V_count,V_count); 
for i = 1:Ni
    idx = newC{i};
    ZC = length(idx);
    for j = 1:(ZC-1)
        adjVV(idx(j),idx(j+1)) = 1;
    end
    adjVV(idx(ZC),idx(1)) = 1;
end
adjVV = adjVV + transpose(adjVV);
adjVV = logical(adjVV);

% under free boundary, the vertex does not always have 3 neighboring cells

% Matrix version of newC
maxzz = max(cellfun('length',newC));
newC_mat = cell2mat(...
                    cellfun(...
                            @(xx)cat(2,xx,zeros(1,maxzz-length(xx))),...
                            newC,'UniformOutput',false...
                            )...
                    );
% Cell-cell adjacency matrix 
adjCC = zeros(Ni,Ni);
for i = 1:V_count
    [m, ~]=find(newC_mat==i);
    if length(m) == 1
        continue;
    elseif length(m) == 2
        adjCC(m(1),m(2))=1;
        adjCC(m(2),m(1))=1;
    elseif length(m) == 3
        adjCC(m(1),m(2))=1;
        adjCC(m(2),m(1))=1;
        adjCC(m(1),m(3))=1;
        adjCC(m(3),m(1))=1;
        adjCC(m(2),m(3))=1;
        adjCC(m(3),m(2))=1;
    end
end

% ------------------------------------------------------------ %
% Vertex-Edge adjacency matrix
% ------------------------------------------------------------ %

[vv1,vv2] = find(triu(adjVV));
nE = length(vv2);

adjVE = zeros(V_count,nE); 
adjVE_oriented = zeros(V_count,nE);
E_list = zeros(nE,2);

for eid = 1:nE
    adjVE(vv1(eid),eid) = 1;
    adjVE(vv2(eid),eid) = 1;
    E_list(eid,1) = vv1(eid);
    E_list(eid,2) = vv2(eid);
    adjVE_oriented(vv1(eid),eid) = -1;
    adjVE_oriented(vv2(eid),eid) = 1;
end

CEadj_oriented = zeros(Ni,nE);
for i = 1:Ni
    vlist = newC{i};
    zz = length(newC{i});
    Vx = newV(vlist,1);
    Vy = newV(vlist,2);
%     [~,ord] = sort(mod(atan2(Vy-mean(Vy),Vx-mean(Vx)),2*pi));
%     vlist = vlist(ord);
    vlist = [vlist, vlist(1)];
    for j = 1:zz
        vid1 = vlist(j);
        vid2 = vlist(j+1);
        eid_candidate = find((E_list(:,1) == vid1) & (E_list(:,2) == vid2));
        if ~isempty(eid_candidate)
            CEadj_oriented(i,eid_candidate) = +1;
        else
            eid_candidate = find((E_list(:,1) == vid2) & (E_list(:,2) == vid1));
            CEadj_oriented(i,eid_candidate) = -1;
        end
    end
end

adjCE = abs(CEadj_oriented);


end