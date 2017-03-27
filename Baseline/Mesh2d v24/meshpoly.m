function [p,t] = meshpoly(node,edge,qtree,p,options)

% MESHPOLY: Core meshing routine called by mesh2d and meshfaces.
%
% Do not call this routine directly, use mesh2d or meshfaces instead!
%
% Inputs:
%
%  NODE     : Nx2 array of geometry XY co-ordinates
%  EDGE     : Mx2 array of connections between NODE, defining geometry 
%             edges
%  QTREE    : Quadtree data structure, defining background mesh and element
%             size function
%  P        : Qx2 array of potential boundary nodes
%  OPTIONS  : Meshing options data structure
%  WBAR     : Handle to progress bar
%
% Outputs:
%
%  P        : Nx2 array of triangle nodes
%  T        : Mx3 array of triangles as indices into P
%
% Mesh2d is a delaunay based algorithm with a "Laplacian-like" smoothing
% operation built into the mesh generation process. 
% 
% An unbalanced quadtree decomposition is used to evaluate the element size 
% distribution required to resolve the geometry. The quadtree is 
% triangulated and used as a backgorund mesh to store the element size 
% data.  
%
% The main method attempts to optimise the node location and mesh topology 
% through an iterative process. In each step a constrained delaunay 
% triangulation is generated with a series of "Laplacian-like" smoothing 
% operations used to improve triangle quality. Nodes are added or removed 
% from the mesh to ensure the required element size distribution is 
% approximated.  
%
% The optimisation process generally returns well shaped meshes with no
% small angles and smooth element size variations. Mesh2d shares some 
% similarities with the Distmesh code: 
%
%   [1] P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
%       SIAM Review, Volume 46 (2), pp. 329-345, June 2004
%
%   Darren Engwirda : 2005-09
%   Email           : d_engwirda@hotmail.com
%   Last updated    : 10/10/2009 with MATLAB 7.0 (Mesh2d v2.4)

shortedge   = 0.75;
longedge    = 1.5;
smalltri    = 0.25;
largetri    = 4.0;
qlimit      = 0.5;
dt          = 0.2;

stats = struct('t_init',0.0,'t_tri',0.0,'t_inpoly',0.0,'t_edge',0.0, ...
                  't_sparse',0.0,'t_search',0.0,'t_smooth',0.0,'t_density',0.0, ...
                     'n_tri',0);

% Initialise mesh
%  P     : Initial nodes
%  T     : Initial triangulation
%  TNDX  : Enclosing triangle for each node as indices into TH
%  FIX   : Indices of FIXED nodes in P
tic
if options.output
   fprintf('Initialising Mesh \n');
end
[p,fix,tndx] = initmesh(p,qtree.p,qtree.t,qtree.h,node,edge);
stats.t_init = toc;

% Main loop
if options.output
   fprintf('Iteration   Convergence (%%)\n');
end
for iter = 1:options.maxit

   [p,i,j] = unique(p,'rows');                                             % Ensure unique node list
   fix = j(fix);
   tndx = tndx(i);
   
   tic
   [p,t] = cdt(p,node,edge);                                               % Constrained Delaunay triangulation
   stats.n_tri = stats.n_tri+1;
   stats.t_tri = stats.t_tri+toc;

   tic
   e = getedges(t,size(p,1));                                              % Unique edges
   stats.t_edge = stats.t_edge+toc;

   % Sparse node-to-edge connectivity matrix
   tic
   nume = size(e,1);
   S = sparse(e(:),[1:nume,1:nume],[ones(nume,1); -ones(nume,1)],size(p,1),nume);
   stats.t_sparse = stats.t_sparse+toc;

   tic
   tndx = mytsearch(qtree.p(:,1),qtree.p(:,2),qtree.t,p(:,1),p(:,2),tndx); % Find enclosing triangle in background mesh for nodes
   hn = tinterp(qtree.p,qtree.t,qtree.h,p,tndx);                           % Size function at nodes via linear interpolation
   h = 0.5*(hn(e(:,1))+hn(e(:,2)));                                        % from the background mesh. Average to edge midpoints.
   stats.t_search = stats.t_search+toc;

   edgev = p(e(:,1),:)-p(e(:,2),:);
   L = max(sqrt(sum((edgev).^2,2)),eps);                                   % Edge lengths 
   
   % Inner smoothing sub-iterations
   time = cputime;
   move = 1.0;
   done = false;
   for subiter = 1:(iter-1)
          
      moveold = move;
    
      % Spring based smoothing
      L0 = h*sqrt(sum(L.^2)/sum(h.^2));
      F = max(L0./L-1.0,-0.1);
      F = S*(edgev.*[F,F]);
      F(fix,:) = 0.0;
      p = p+dt*F;

      % Measure convergence
      edgev = p(e(:,1),:)-p(e(:,2),:);
      L0 = max(sqrt(sum((edgev).^2,2)),eps);                               % Edge lengths 
      move = norm((L0-L)./L,'inf');                                        % Percentage change
      L = L0;

      if move<options.mlim                                                 % Test convergence
         done = true;
         break         
      end
      
   end
   stats.t_smooth = stats.t_smooth+(cputime-time);
   
   if options.output
      fprintf('%2i           %2.1f\n',[iter, 100.0*min(1.0, options.mlim/max(move,eps))]);
   end
   
   tic
   [p,t] = cdt(p,node,edge);                                               % Constrained Delaunay triangulation
   stats.n_tri = stats.n_tri+1;
   stats.t_tri = stats.t_tri+toc;

   tic
   e = getedges(t,size(p,1));                                              % Unique edges
   stats.t_edge = stats.t_edge+toc;

   edgev = p(e(:,1),:)-p(e(:,2),:);
   L = max(sqrt(sum((edgev).^2,2)),eps);                                   % Edge lengths 

   tic
   tndx = mytsearch(qtree.p(:,1),qtree.p(:,2),qtree.t,p(:,1),p(:,2),tndx); % Find enclosing triangle in background mesh for nodes
   hn = tinterp(qtree.p,qtree.t,qtree.h,p,tndx);                           % Size function at nodes via linear interpolation
   h = 0.5*(hn(e(:,1))+hn(e(:,2)));                                        % from the background mesh. Average to edge midpoints.
   stats.t_search = stats.t_search+toc;
   
   r = L./h;
   if done && (max(r)<3.0)                                                 % Main loop convergence
      break
   else
      if (iter==options.maxit)
         fprintf('WARNING: Maximum number of iterations reached. Solution did not converge! \n');
         fprintf('Please email the geometry and settings to d_engwirda@hotmail.com \n');
      end
   end
   
   % Nodal density control
   tic
   if iter<options.maxit      
      % Estimate required triangle area from size function
      Ah = 0.5*tricentre(t,hn).^2;
      % Remove nodes
      i = find(abs(triarea(p,t))<smalltri*Ah);                             % Remove all nodes in triangles with small area
      k = find(sum(abs(S),2)<2);                                           % Nodes with less than 2 edge connections
      j = find(r<shortedge);                                               % Remove both nodes for short edges
      if ~isempty(j) || ~isempty(k) || ~isempty(i)
         prob = false(size(p,1),1);                                        % True for nodes to be removed
         prob(e(j,:)) = true;                                              % Edges with insufficient length
         prob(t(i,:)) = true;                                              % Triangles with insufficient area
         prob(k) = true;                                                   % Remove nodes with less than 2 edge connections
         prob(fix) = false;                                                % Don't remove fixed nodes
         pnew = p(~prob,:);                                                % New node list
         tndx = tndx(~prob);        
         j = zeros(size(p,1),1);                                           % Re-index FIX to keep consistent
         j(~prob) = 1;
         j = cumsum(j);
         fix = j(fix);
      else
         pnew = p;                                                  
      end
      % Add new nodes 
      i = abs(triarea(p,t))>largetri*Ah;                                   % Large triangles
      r = longest(p,t)./tricentre(t,hn);
      k = (r>longedge) & (quality(p,t)<qlimit);                            % Low quality triangles
      if any(i|k)

         k = find(k & ~i);
         i = find(i);
         
         % Add new nodes at circumcentres
         cc = circumcircle(p,[t(i,:); t(k,:)]);
         
         % Don't add multiple points in one circumcircle
         ok = [true(size(i)); false(size(k))];
         for ii = (length(i)+1):size(cc,1)
            % Current point
            x = cc(ii,1);
            y = cc(ii,2);
            % Check if inside any accepted circumcircles
            in = false;
            j = find(ok);
            for jj = 1:length(j)
               kk = j(jj);
               dx = (x-cc(kk,1))^2;
               if dx<cc(kk,3) && (dx+(y-cc(kk,2))^2)<cc(kk,3)
                  in = true;
                  break;
               end
            end
            if ~in
               ok(ii) = true;
            end
         end
         cc = cc(ok,:);
         cc = cc(inpoly(cc(:,1:2),node,edge),:);                           % Only take internal points
         
         % New arrays
         pnew = [pnew; cc(:,1:2)];
         tndx = [tndx; zeros(size(cc,1),1)];                               
      end
      p = pnew;
   end
   stats.t_density = stats.t_density+toc;

end

if options.debug
   disp(stats);
end

end      % meshpoly()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,t] = cdt(p,node,edge)

% Approximate geometry-constrained Delaunay triangulation.

t = mydelaunayn(p);                                                        % Delaunay triangulation via QHULL

% Impose geometry constraints
i = inpoly(tricentre(t,p),node,edge);                                      % Take triangles with internal centroids
t = t(i,:);

end      % [p,t] = cdt(p,node,edge)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,fix,tndx] = initmesh(p,ph,th,hh,node,edge)

% Initialise the mesh nodes

% Boundary nodes for all geometry edges have been passed in. Only take
% those in the current face
i = findedge(p,node,edge,1.0e-08);
p = p(i>0,:);
fix = (1:size(p,1))';

% Initial nodes taken as fixed boundary nodes + internal nodes from
% quadtree.
[i,j] = inpoly(ph,node,edge);
p = [p; ph(i&~j,:)];
tndx = zeros(size(p,1),1);

end      % initmesh()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = getedges(t,n)

% Get the unique edges and boundary nodes in a triangulation.

e = sortrows( sort([t(:,[1,2]); t(:,[1,3]); t(:,[2,3])],2) );
idx = all(diff(e,1)==0,2);                                                 % Find shared edges
idx = [idx;false]|[false;idx];                                             % True for all shared edges
bnd = e(~idx,:);                                                           % Boundary edges
e = e(idx,:);                                                              % Internal edges
e = [bnd; e(1:2:end-1,:)];                                                 % Unique edges and bnd edges for tri's

end      % getedges()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc = tricentre(t,f)

% Interpolate nodal F to the centroid of the triangles T.

fc = (f(t(:,1),:)+f(t(:,2),:)+f(t(:,3),:))/3.0;

end      % tricentre()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = longest(p,t)

% Return the length of the longest edge in each triangle.

d1 = sum((p(t(:,2),:)-p(t(:,1),:)).^2,2);
d2 = sum((p(t(:,3),:)-p(t(:,2),:)).^2,2);
d3 = sum((p(t(:,1),:)-p(t(:,3),:)).^2,2);

d = sqrt(max([d1,d2,d3],[],2));

end      % longest()
