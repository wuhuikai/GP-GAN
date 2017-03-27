function fi = tinterp(p,t,f,pi,i)

%  TINTERP: Triangle based linear interpolation.
%
%  fi = tinterp(p,t,f,pi,i);
%
%  p  : Nx2 array of nodal XY coordinates, [x1,y1; x2,y2; etc]
%  t  : Mx3 array of triangles as indices, [n11,n12,n13; n21,n22,n23; etc]
%  f  : Nx1 function vector, f(x,y)
%  pi : Jx2 matrix of interpolation points
%  fi : Jx1 interpolant function vector, fi(xi,yi)
%
% Performs nearest-neighbour extrapolation for points outside the
% triangulation.

% Darren Engwirda - 2005-2007

% Alloc output
fi = zeros(size(pi,1),1);

% Deal with points oustide convex hull
out = isnan(i);
if any(out)
   % Do nearest neighbour extrapolation for outside points
   d = dsearchn(p(:,1),p(:,2),t,pi(out,1),pi(out,2));
   fi(out) = f(d);
end




% Keep internal points
pin = pi(~out,:);
tin = t(i(~out),:);

% Corner nodes
t1 = tin(:,1);
t2 = tin(:,2);
t3 = tin(:,3);

% Calculate areas
dp1 = pin-p(t1,:);
dp2 = pin-p(t2,:);
dp3 = pin-p(t3,:);
A3 = abs(dp1(:,1).*dp2(:,2)-dp1(:,2).*dp2(:,1));
A2 = abs(dp1(:,1).*dp3(:,2)-dp1(:,2).*dp3(:,1));
A1 = abs(dp3(:,1).*dp2(:,2)-dp3(:,2).*dp2(:,1));

% Linear interpolation
fi(~out) = (A1.*f(t1)+A2.*f(t2)+A3.*f(t3))./(A1+A2+A3);

end      % tinterp()
