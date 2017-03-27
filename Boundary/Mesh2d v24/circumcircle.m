function cc = circumcircle(p,t)

% CIRCUMCIRCLE: XY centre co-ordinates and radius of triangle 
% circumcircles.
%
% P   : Nx2 array of nodal XY co-ordinates
% T   : Mx3 array of triangles as indices into P
% CC  : Mx3 array of circimcircles CC(:,1:2) = XY, CC(:,3) = R^2

cc = zeros(size(t));

% Corner XY
p1 = p(t(:,1),:); 
p2 = p(t(:,2),:); 
p3 = p(t(:,3),:);

% Set equation for center of each circumcircle: 
%    [a11,a12; a21,a22] * [x; y] = [b1; b2] * 0.5;
a1 = p2-p1; 
a2 = p3-p1; 
b1 = sum(a1.*(p2+p1),2); %a1(:,1).*(p2(:,1)+p1(:,1)) + a1(:,2).*(p2(:,2)+p1(:,2));
b2 = sum(a2.*(p3+p1),2); %a2(:,1).*(p3(:,1)+p1(:,1)) + a2(:,2).*(p3(:,2)+p1(:,2));

% Explicit inversion
idet = 0.5./(a1(:,1).*a2(:,2)-a2(:,1).*a1(:,2));

% Circumcentre XY
cc(:,1) = ( a2(:,2).*b1 - a1(:,2).*b2).*idet;
cc(:,2) = (-a2(:,1).*b1 + a1(:,1).*b2).*idet;

% Radius^2
cc(:,3) = sum((p1-cc(:,1:2)).^2,2);

end      % circumcircle()