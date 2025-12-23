function [HxAx, HxAy, HxBx, HxBy, HxCx, HxCy,...
          HyAx, HyAy, HyBx, HyBy, HyCx, HyCy] = get_dH(Ax,Ay,Bx,By,Cx,Cy)

% partial derivative of  vertex position with respect to cell center positions;

HxAx = (2*Ax*(By - Cy))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy)) - ((2*By - 2*Cy)*((By - Cy)*(Ax^2 + Ay^2) - (Ay - Cy)*(Bx^2 + By^2) + (Ay - By)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2;
HxAy = (Cx^2 - By^2 - Bx^2 + Cy^2 + 2*Ay*(By - Cy))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy)) + ((2*Bx - 2*Cx)*((By - Cy)*(Ax^2 + Ay^2) - (Ay - Cy)*(Bx^2 + By^2) + (Ay - By)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2;
HxBx = ((2*Ay - 2*Cy)*((By - Cy)*(Ax^2 + Ay^2) - (Ay - Cy)*(Bx^2 + By^2) + (Ay - By)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2 - (2*Bx*(Ay - Cy))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy));
HxBy = - (Cx^2 - Ay^2 - Ax^2 + Cy^2 + 2*By*(Ay - Cy))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy)) - ((2*Ax - 2*Cx)*((By - Cy)*(Ax^2 + Ay^2) - (Ay - Cy)*(Bx^2 + By^2) + (Ay - By)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2;
HxCx = (2*Cx*(Ay - By))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy)) - ((2*Ay - 2*By)*((By - Cy)*(Ax^2 + Ay^2) - (Ay - Cy)*(Bx^2 + By^2) + (Ay - By)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2;
HxCy = (Bx^2 - Ay^2 - Ax^2 + By^2 + 2*Cy*(Ay - By))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy)) + ((2*Ax - 2*Bx)*((By - Cy)*(Ax^2 + Ay^2) - (Ay - Cy)*(Bx^2 + By^2) + (Ay - By)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2;
HyAx = ((2*By - 2*Cy)*((Bx - Cx)*(Ax^2 + Ay^2) - (Ax - Cx)*(Bx^2 + By^2) + (Ax - Bx)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2 - (Cx^2 - By^2 - Bx^2 + Cy^2 + 2*Ax*(Bx - Cx))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy));
HyAy = - (2*Ay*(Bx - Cx))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy)) - ((2*Bx - 2*Cx)*((Bx - Cx)*(Ax^2 + Ay^2) - (Ax - Cx)*(Bx^2 + By^2) + (Ax - Bx)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2;
HyBx = (Cx^2 - Ay^2 - Ax^2 + Cy^2 + 2*Bx*(Ax - Cx))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy)) - ((2*Ay - 2*Cy)*((Bx - Cx)*(Ax^2 + Ay^2) - (Ax - Cx)*(Bx^2 + By^2) + (Ax - Bx)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2;
HyBy = (2*By*(Ax - Cx))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy)) + ((2*Ax - 2*Cx)*((Bx - Cx)*(Ax^2 + Ay^2) - (Ax - Cx)*(Bx^2 + By^2) + (Ax - Bx)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2;
HyCx = ((2*Ay - 2*By)*((Bx - Cx)*(Ax^2 + Ay^2) - (Ax - Cx)*(Bx^2 + By^2) + (Ax - Bx)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2 - (Bx^2 - Ay^2 - Ax^2 + By^2 + 2*Cx*(Ax - Bx))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy));
HyCy = - (2*Cy*(Ax - Bx))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy)) - ((2*Ax - 2*Bx)*((Bx - Cx)*(Ax^2 + Ay^2) - (Ax - Cx)*(Bx^2 + By^2) + (Ax - Bx)*(Cx^2 + Cy^2)))/(2*Cx*(Ay - By) - 2*Bx*(Ay - Cy) + 2*Ax*(By - Cy))^2;






end