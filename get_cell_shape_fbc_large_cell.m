% Free Boundary Conditions version
function [perim_list, area_list] = get_cell_shape_fbc_large_cell(A_total,V,newC,cell_large)

        Ni = size(newC,1);
        Npts = Ni + 1;
        area_list = zeros(Npts,1);
        perim_list = zeros(Npts,1);
        for i = 1:Ni
            Ximg = V(newC{i},:);
            Z = size(Ximg,1);
%             [~,ord] = sort(mod(atan2(Ximg(:,2)-mean(Ximg(:,2)),Ximg(:,1)-mean(Ximg(:,1))),2*pi));
%             Ximg = Ximg(ord,:);
            dx = Ximg( [ 2:Z 1 ], 1) - Ximg(:,1);
            dy = Ximg( [ 2:Z 1 ], 2) - Ximg(:,2);

            % summations for CW boundary integrals
            area_list(i) = abs(sum( Ximg(:,2).*dx - Ximg(:,1).*dy )/2);
            perim_list(i) = sum( sqrt( dx.*dx +dy.*dy ) );
        end

        % area of the large cell;
        area_list(Npts) = (A_total - sum(area_list(1:Ni)));

%         if sum(area_list(1:Ni)) > A_total
%             area_list(Npts) = 1;
%         end

        % perimeter of the large cell;
        Ximg = V(cell_large,:);
        Z = size(Ximg,1);
%         [~,ord] = sort(mod(atan2(Ximg(:,2)-mean(Ximg(:,2)),Ximg(:,1)-mean(Ximg(:,1))),2*pi));
%         Ximg = Ximg(ord,:);
        dx = Ximg( [ 2:Z 1 ], 1) - Ximg(:,1);
        dy = Ximg( [ 2:Z 1 ], 2) - Ximg(:,2);
        perim_list(Npts) = sum( sqrt( dx.*dx +dy.*dy ) );
end
