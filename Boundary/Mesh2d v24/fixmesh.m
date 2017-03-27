function [p,t,pfun,tfun] = fixmesh(p,t,pfun,tfun)

%  FIXMESH: Ensure that triangular mesh data is consistent.
%
%  [p,t,pfun,tfun] = fixmesh(p,t,pfun,tfun);
%
%  p     : Nx2 array of nodal XY coordinates, [x1,y1; x2,y2; etc]
%  t     : Mx3 array of triangles as indices, [n11,n12,n13; n21,n22,n23;
%          etc]
%  pfun  : (Optional) NxK array of nodal function values. Each column in
%          PFUN corresponds to a dependent function, PFUN(:,1) = F1(P),
%          PFUN(:,2) = F2(P) etc, defined at the nodes.
%  tfun  : (Optional) MxK array of triangle function values. Each column in
%          TFUN corresponds to a dependent function, TFUN(:,1) = F1(T),
%          TFUN(:,2) = F2(T) etc, defined on the triangles.
%
% The following checks are performed:
%
%  1. Nodes not refereneced in T are removed.
%  2. Duplicate nodes are removed.
%  3. Triangles are ordered counter-clockwise.
%  4. Triangles with an area less than 1.0e-10*eps*norm(A,'inf')
%     are removed

% Darren Engwirda - 2007.

TOL = 1.0e-10;

if (nargin<4)
   tfun = [];
   if (nargin<3)
      pfun = [];
      if nargin<2
         error('Wrong number of inputs');
      end
   end
elseif (nargin>4)
   error('Wrong number of inputs');
end
if (nargout>4)
   error('Wrong number of outputs');
end
if (numel(p)~=2*size(p,1))
   error('P must be an Nx2 array');
end
if (numel(t)~=3*size(t,1))
   error('T must be an Mx3 array');
end
if (any(t(:))<1) || (max(t(:))>size(p,1))
   error('Invalid T');
end
if ~isempty(pfun)
   if (size(pfun,1)~=size(p,1)) || (ndims(pfun)~=2)
      error('PFUN must be an NxK array');
   end
end
if ~isempty(tfun)
   if (size(tfun,1)~=size(t,1)) || (ndims(tfun)~=2)
      error('TFUN must be an Mxk array');
   end
end

% Remove duplicate nodes
[i,i,j] = unique(p,'rows');
if ~isempty(pfun)
   pfun = pfun(i,:);
end
p = p(i,:);
t = reshape(j(t),size(t));

% Triangle area
A = triarea(p,t);
Ai = A<0.0;
Aj = abs(A)>TOL*norm(A,'inf');

% Flip node numbering to give a counter-clockwise order
t(Ai,[1,2]) = t(Ai,[2,1]);

% Remove zero area triangles
t = t(Aj,:);
if ~isempty(tfun)
   tfun = tfun(Aj,:);
end

% Remove un-used nodes
[i,j,j] = unique(t(:));
if ~isempty(pfun)
   pfun = pfun(i,:);
end
p = p(i,:);
t = reshape(j,size(t));

end      % fixmesh()
