function A = triarea(p,t)

% TRIAREA: Area of triangles assuming counter-clockwise (CCW) node
% ordering.
%
%  P  : Nx2 array of XY node co-ordinates
%  T  : Mx3 array of triangles as indices into P
%  A  : Mx1 array of triangle areas

% Darren Engwirda - 2007

d12 = p(t(:,2),:)-p(t(:,1),:);
d13 = p(t(:,3),:)-p(t(:,1),:);
A = (d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1));

end      % triarea()