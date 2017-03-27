function [p,t] = smoothmesh(p,t,maxit,tol)

%  SMOOTHMESH: Smooth a triangular mesh using Laplacian smoothing.
%
% Laplacian smoothing is an iterative process that generally leads to an
% improvement in the quality of the elements in a triangular mesh.
%
%  [p,t] = smoothmesh(p,t);
%
%  p     : Nx2 array of nodal XY coordinates, [x1,y1; x2,y2; etc].
%  t     : Mx3 array of triangles as indices, [n11,n12,n13; 
%                                              n21,n22,n23; etc].
%  maxit : Maximum allowable iterations.
%  tol   : Convergence tolerance (Percentage change in edge length must be 
%          less than TOL).
%
% If MAXIT or TOL are left empty the default values MAXIT = 20 and TOL =
% 0.01 are used.
%
%  EXAMPLE:
%
%  [p,t] = smoothmesh(p,t,10,0.05);
%
% See also MESH2D, REFINE

% Darren Engwirda - 2007.

if nargin<4
   tol = [];
   if nargin<3
      maxit = [];
      if nargin<2
         error('Incorrect number of inputs.');
      end
   end
elseif nargin>5
   error('Incorrect number of inputs.')
end
if nargout>2
   error('Incorrect number of outputs.');
end
if isempty(tol)
   tol = 0.01;
end
if isempty(maxit)
   maxit = 20;
end

[p,t] = fixmesh(p,t);                                                      % Ensure consistent mesh

n = size(p,1);
S = sparse(t(:,[1,1,2,2,3,3]),t(:,[2,3,1,3,1,2]),1,n,n);                   % Sparse connectiity matrix
W = sum(S,2);
if any(W==0)
   error('Invalid mesh. Hanging nodes found.');
end

edge = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])];                               % Non-unique edges
edge = sortrows(sort(edge,2));                                             % Put shared edges next to each other
idx = all(diff(edge)==0,2);                                                % Find shared edges
idx = [idx;false]|[false;idx];                                             % True for all shared edges
bnde = edge(~idx,:);                                                       % Boundary edges
edge = edge(idx,:);                                                        % Internal edges
edge = [bnde; edge(1:2:end-1,:)];                                          % Unique edges
bnd = unique(bnde);                                                        % Boundary nodes

L = max(sqrt(sum((p(edge(:,1),:)-p(edge(:,2),:)).^2,2)),eps);              % Edge length

for iter = 1:maxit
   pnew = (S*p)./[W,W];                                                    % Laplacian smoothing
   pnew(bnd,:) = p(bnd,:);                                                 % Dont let BND nodes move
   p = pnew;

   Lnew = max(sqrt(sum((p(edge(:,1),:)-p(edge(:,2),:)).^2,2)),eps);        % Edge length
   move = norm((Lnew-L)./Lnew,inf);                                        % Percentage change
   if move<tol
      break
   end
   L = Lnew;
end
if iter==maxit
   disp('WARNING: Maximum number of iterations reached, solution did not converge!');
end

end      % smoothmesh()
