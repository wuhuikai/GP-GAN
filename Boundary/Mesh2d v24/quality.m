function q = quality(p,t)

%  QUALITY: Approximate triangle quality. 
%
%  q = quality(p,t);
%
%  p: Nx2 array of nodal XY coordinates, [x1,y1; x2,y2; etc]
%  t: Mx3 array of triangles as indices, [n11,n12,n13; n21,n22,n23; etc]
%  q: Mx1 vector of triangle qualities. 0<=q<=1.

% Darren Engwirda - 2007.

% Nodes
p1 = p(t(:,1),:); 
p2 = p(t(:,2),:); 
p3 = p(t(:,3),:);

% Approximate quality
d12 = p2-p1;
d13 = p3-p1;
d23 = p3-p2;
q = 3.4641*abs(d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))./sum(d12.^2+d13.^2+d23.^2,2);

end      % quality()
