function i = mytsearch(x,y,t,xi,yi,i)

%  MYTSEARCH: Find the enclosing triangle for points in a 2D plane.
%
%  i = mytsearch(x,y,t,xi,yi,iguess);
%
% The indices of the triangles enclosing the points in [XI,YI] are
% returned. The triangulation T of [X,Y] must be convex. Points lying
% outside the triangulation are assigned a NaN index.
%
% IGUESS is an optional initial guess for the indicies. A full search is
% done using the standard TSEARCH function for points with an invalid
% initial guess.

% Darren Engwirda - 2007.

% I/O and error checks
if nargin<6
   i = [];
   if nargin<5
      error('Wrong number of inputs');
   end
elseif nargin>6
   error('Wrong number of inputs');
end
if nargout>1
   error('Wrong number of outputs');
end
ni = size(xi,1);
if (numel(xi)~=ni) || (numel(yi)~=ni) || ...
      (numel(x)~=numel(y)) || (~isempty(i) && (numel(i)~=ni))
   error('Wrong input dimensions');
end

% Translate to the origin and scale the min xy range onto [-1,1]
% This is absolutely critical to avoid precision issues for large problems!
maxxy = max([x,y]);
minxy = min([x,y]);
den = 0.5*min(maxxy-minxy);

x  = ( x-0.5*(minxy(1)+maxxy(1))) / den;
y  = ( y-0.5*(minxy(2)+maxxy(2))) / den;
xi = (xi-0.5*(minxy(1)+maxxy(1))) / den;
yi = (yi-0.5*(minxy(2)+maxxy(2))) / den;

% Check initial guess
if ~isempty(i)
   k = find(i>0 & ~isnan(i));

   tri = i(k);

   n1 = t(tri,1);
   n2 = t(tri,2);
   n3 = t(tri,3);

   ok = sameside(x(n1),y(n1),x(n2),y(n2),xi(k),yi(k),x(n3),y(n3)) & ...
        sameside(x(n2),y(n2),x(n3),y(n3),xi(k),yi(k),x(n1),y(n1)) & ...
        sameside(x(n3),y(n3),x(n1),y(n1),xi(k),yi(k),x(n2),y(n2));

   j = true(ni,1);
   j(k(ok)) = false;
else
   j = true(ni,1);
end

% Do a full search for points that failed
if any(j)
   i(j) = tsearchn([x,y],t,[xi(j),yi(j)]);
end

end      % mytsearch()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function i = sameside(xa,ya,xb,yb,x1,y1,x2,y2)

% Test if [x1(i),y1(i)] and [x2(i),y2(i)] lie on the same side of the line
% AB(i).

dx = xb-xa;
dy = yb-ya;
a1 = (x1-xa).*dy-(y1-ya).*dx;
a2 = (x2-xa).*dy-(y2-ya).*dx;

% If sign(a1)=sign(a2) the points lie on the same side
i = false(length(xa),1);
i(a1.*a2>=0.0) = true;

end      % sameside()
