function [e,te,e2t,bnd] = connectivity(p,t)

%  CONNECTIVITY: Assemble connectivity data for a triangular mesh.
%
% The edge based connectivity is built for a triangular mesh and the
% boundary nodes identified. This data should be useful when implementing
% FE/FV methods using triangular meshes.
%
%  [e,te,et2,bnd] = connectivity(p,t);
%
%  p   : Nx2 array of nodes coordinates, [x1,y1; x2,y2; etc]
%  t   : Mx3 array of triangles as indices, [n11,n12,n13; n21,n22,n23; etc]
%
%  e   : Kx2 array of unique mesh edges - [n11,n12; n21,n22; etc]
%  te  : Mx3 array of triangles as indices into E, [e11,e12,e13; 
%                                                   e21,e22,e23; etc]
%  e2t : Kx2 array of triangle neighbours for unique mesh edges -
%        [t11,t12; t21,t22; etc]. Each row has two entries corresponding to
%        the triangle numbers associated with each edge in E. Boundary
%        edges have e2t(i,2)=0.
%  bnd : Nx1 logical array identifying boundary nodes. P(i,:) is a boundary
%        node if BND(i)=TRUE.
%
% See also MESH2D, REFINE

% Darren Engwirda - 2007

if nargin<2
   error('Wrong number of inputs');
end
if nargout>4
   error('Wrong number of outputs');
end
if numel(p)~=2*size(p,1)
   error('P must be an Nx2 array');
end
if numel(t)~=3*size(t,1)
   error('T must be an Mx3 array');
end
if any(t(:)<1) || max(t(:)>size(p,1))
   error('Invalid T');
end

% Unique mesh edges as indices into P
numt = size(t,1);
vect = 1:numt;                                                             % Triangle indices
e = [t(:,[1,2]); t(:,[2,3]); t(:,[3,1])];                                  % Edges - not unique
[e,j,j] = unique(sort(e,2),'rows');                                        % Unique edges
te = [j(vect), j(vect+numt), j(vect+2*numt)];                              % Unique edges in each triangle

% Edge-to-triangle connectivity
% Each row has two entries corresponding to the triangle numbers
% associated with each edge. Boundary edges have e2t(i,2)=0.
nume = size(e,1);
e2t  = zeros(nume,2);
for k = 1:numt
   j = 1;
   while j<=3
      ce = te(k,j);
      if e2t(ce,1)==0
         e2t(ce,1) = k;
      else
         e2t(ce,2) = k;
      end
      j = j+1;
   end
end

% Flag boundary nodes
bnd = false(size(p,1),1);
bnd(e(e2t(:,2)==0,:)) = true;                                              % True for bnd nodes

end      % connectivity()
