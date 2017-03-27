function [p,t,h] = quadtree(node,edge,hdata,dhmax,output)

%  QUADTREE: 2D quadtree decomposition of polygonal geometry.
%
% The polygon is first rotated so that the minimal enclosing rectangle is
% aligned with the Cartesian XY axes. The long axis is aligned with Y. This
% ensures that the quadtree generated for a geometry input that has
% undergone arbitrary rotations in the XY plane is always the same.
%
% The bounding box is recursively subdivided until the dimension of each 
% box matches the local geometry feature size. The geometry feature size is 
% based on the minimum distance between linear geometry segments.
%
% A size function is obtained at the quadtree vertices based on the minimum
% neighbouring box dimension at each vertex. This size function is gradient
% limited to produce a smooth function.
%
% The quadtree is triangulated to form a background mesh, such that the
% size function may be obtained at any XY position within the domain via
% triangle based linear interpolation. The triangulation is done based on
% the quadtree data structures directly (i.e. NOT using Qhull) which is 
% more reliable and produces a consistently oriented triangulation.
%
% The initial rotations are undone.
%
%  node  : [x1,y1; x2,y2; etc] geometry vertices
%  edge  : [n11,n12; n21,n22; etc] geometry edges as connections in NODE
%  hdata : User defined size function structure
%  dhmax : Maximum allowalble relative gradient in the size function
%  wbar  : Handle to progress bar from MESH2D
%
%   p    : Background mesh nodes
%   t    : Background mesh triangles
%   h    : Size function value at p

%   Darren Engwirda : 2007
%   Email           : d_engwirda@hotmail.com
%   Last updated    : 18/11/2007 with MATLAB 7.0

% Bounding box
XYmax = max(node,[],1);
XYmin = min(node,[],1);

% Rotate NODE so that the long axis of the minimum bounding rectangle is
% aligned with the Y axis.
theta = minrectangle(node);
node = rotate(node,theta);

% Rotated XY edge endpoints
edgexy = [node(edge(:,1),:), node(edge(:,2),:)];

% LOCAL FEATURE SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if output
   fprintf('Estimating local geometry feature size\n');
end

% Get size function data
[hmax,edgeh,fun,args] = gethdata(hdata);

% Insert test points along the boundaries at which the LFS can be
% approximated.
wm = 0.5*(edgexy(:,[1,2])+edgexy(:,[3,4]));                                % Use the edge midpoints as a first pass
len = sqrt(sum((edgexy(:,[3,4])-edgexy(:,[1,2])).^2,2));                   % Edge length
L = 2.0*dist2poly(wm,edgexy,2.0*len);                                      % Estimate the LFS at these points by calculating
                                                                           % the distance to the closest edge segment                           
% In cases where edges are separated by less than their length
% we will need to add more points to capture the LFS in these
% regions. This allows us to pick up "long and skinny" geometry
% features
r = 2.0*len./L;                                                            % Compare L (LFS approximation at wm) to the edge lengths
r = round((r-2.0)/2.0);                                                    % Amount of points that need to be added
add = find(r);                                                             % at each edge
if ~isempty(add)
   num = 2*sum(r(add));                                                    % Total number of points to be added
   start = size(wm,1)+1;
   wm = [wm; zeros(num,2)];                                                % Alloc space
   L = [L; zeros(num,1)];
   next = start;
   for j = 1:length(add)                                                   % Loop through edges to be subdivided
      
      ce = add(j);                                                         % Current edge
      num = r(ce);
      tmp = (1:num)'/(num+1);                                              % Subdivision increments
      num = next+2*num-1;

      x1 = edgexy(ce,1); x2 = edgexy(ce,3); xm = wm(ce,1);                 % Edge values
      y1 = edgexy(ce,2); y2 = edgexy(ce,4); ym = wm(ce,2);

      wm(next:num,:) = [x1+tmp*(xm-x1), y1+tmp*(ym-y1)                     % Add to list
                        xm+tmp*(x2-xm), ym+tmp*(y2-ym)];
      
      L(next:num) = L(ce);                                                 % Upper bound on LFS
                     
      next = num+1;
   
   end
   L(start:end) = dist2poly(wm(start:end,:),edgexy,L(start:end));          % Estimate LFS at the new points
end

% Compute the required size along the edges for any boundary layer size
% functions and add additional points where necessary.
if ~isempty(edgeh)
   for j = 1:size(edgeh,1)
      if L(edgeh(j,1))>=edgeh(j,2)
        
         cw = edgeh(j,1);
         r = 2.0*len(cw)/edgeh(j,2);
         r = ceil((r)/2.0);                                          % Number of points to be added
         tmp = (1:r-1)'/r;

         x1 = edgexy(cw,1); x2 = edgexy(cw,3); xm = wm(cw,1);              % Edge values
         y1 = edgexy(cw,2); y2 = edgexy(cw,4); ym = wm(cw,2);

         wm = [wm; x1+tmp*(xm-x1), y1+tmp*(ym-y1);                         % Add to list
                   xm+tmp*(x2-xm), ym+tmp*(y2-ym)];

         L(cw) = edgeh(j,2);                                                    % Update LFS
         L = [L; edgeh(j,2)*ones(2*length(tmp),1)];
      
      end
   end
end

% To speed the point location in the quadtree decomposition
% sort the LFS points based on y-value
[i,i] = sort(wm(:,2));
wm = wm(i,:);
L = L(i);
nw = size(wm,1);

% UNBALANCED QUADTREE DECOMPOSITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if output
   fprintf('Quadtree decomposition\n');
end

xymin = min([edgexy(:,[1,2]); edgexy(:,[3,4])]);                           % Bounding box
xymax = max([edgexy(:,[1,2]); edgexy(:,[3,4])]);

dim = 2.0*max(xymax-xymin);                                                    % Bbox dimensions
xm = 0.5*(xymin(1)+xymax(1));
ym = 0.5*(xymin(2)+xymax(2));

% Setup boxes with a consistent CCW node order
%  b(k,1) = n1 : bottom left
%  b(k,2) = n2 : bottom right
%  b(k,3) = n3 : top right
%  b(k,4) = n4 : top left

% Start with bbox
p = [xm-0.5*dim, ym-0.5*dim
     xm+0.5*dim, ym-0.5*dim
     xm+0.5*dim, ym+0.5*dim
     xm-0.5*dim, ym+0.5*dim];
b = [1,2,3,4];

% User defined size functions
pr = rotate(p,-theta);
h = userhfun(pr(:,1),pr(:,2),fun,args,hmax,XYmin,XYmax);

pblock = 5*nw;                                                             % Alloc memory in blocks
bblock = pblock;

np = size(p,1);
nb = size(b,1);
test = true(nb,1);
while true                                                           
   
   vec = find(test(1:nb));                                                 % Boxes to be checked at this step
   if isempty(vec)
      break
   end

   N = np;
   for k = 1:length(vec)                                                   % Loop through boxes to be checked for subdivision
      
      m  = vec(k);                                                         % Current box

      n1 = b(m,1);   n2 = b(m,2);                                          % Corner nodes
      n3 = b(m,3);   n4 = b(m,4);
      x1 = p(n1,1);  y1 = p(n1,2);                                         % Corner xy
      x2 = p(n2,1);  y4 = p(n4,2);

      % Binary search to find first wm with y>=ymin for current box
      if wm(1,2)>=y1
         start = 1;
      elseif wm(nw,2)<y1
         start = nw+1;
      else
         lower = 1;
         upper = nw;
         for i = 1:nw
            start = round(0.5*(lower+upper));
            if wm(start,2)<y1
               lower = start;
            elseif wm(start-1,2)<y1
               break;
            else
               upper = start;
            end
         end
      end
      
      % Init LFS as the min of corner user defined size function values
      LFS = 1.5*h(n1);
      if 1.5*h(n2)<LFS, LFS = 1.5*h(n2); end
      if 1.5*h(n3)<LFS, LFS = 1.5*h(n3); end
      if 1.5*h(n4)<LFS, LFS = 1.5*h(n4); end

      % Loop through all WM in box and take min LFS
      for i = start:nw                                                     % Loop through points (acending y-value order)
         if (wm(i,2)<=y4)                                                  % Check box bounds and current min
            if (wm(i,1)>=x1) && (wm(i,1)<=x2) && (L(i)<LFS)
               LFS = L(i);                                                 % New min found - reset
            end
         else                                                              % Due to the sorting
            break;
         end
      end

      % Split box into 4
      if (x2-x1)>=LFS 
         
         if (np+5)>=size(p,1)                                              % Alloc memory on demand
            p = [p; zeros(pblock,2)];
            pblock = 2*pblock;
         end
         if (nb+3)>=size(b,1)
            b = [b; zeros(bblock,4)];
            test = [test; true(bblock,1)];
            bblock = 2*bblock;
         end

         xm = x1+0.5*(x2-x1);                                              % Current midpoints
         ym = y1+0.5*(y4-y1);

         % New nodes
         p(np+1,1) = xm;   p(np+1,2) = ym;
         p(np+2,1) = xm;   p(np+2,2) = y1;
         p(np+3,1) = x2;   p(np+3,2) = ym;
         p(np+4,1) = xm;   p(np+4,2) = y4;
         p(np+5,1) = x1;   p(np+5,2) = ym;   

         % New boxes
         b(m,1)      = n1;               % Box 1
         b(m,2)      = np+2;
         b(m,3)      = np+1;
         b(m,4)      = np+5;
         b(nb+1,1)   = np+2;             % Box 2
         b(nb+1,2)   = n2;
         b(nb+1,3)   = np+3;
         b(nb+1,4)   = np+1;
         b(nb+2,1)   = np+1;             % Box 3
         b(nb+2,2)   = np+3;
         b(nb+2,3)   = n3;
         b(nb+2,4)   = np+4;
         b(nb+3,1)   = np+5;             % Box 4
         b(nb+3,2)   = np+1;
         b(nb+3,3)   = np+4;
         b(nb+3,4)   = n4;

         nb = nb+3;
         np = np+5;
      else
         test(m) = false;
      end
   end

   % User defined size function at new nodes
   pr = rotate(p(N+1:np,:),-theta);
   h = [h; userhfun(pr(:,1),pr(:,2),fun,args,hmax,XYmin,XYmax)];

end
p = p(1:np,:);
b = b(1:nb,:);

% Remove duplicate nodes
[p,i,j] = unique(p,'rows');                               
h = h(i);
b = reshape(j(b),size(b));

% FORM SIZE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if output
   fprintf('Forming element size function\n');
end

% Unique edges
e = unique(sort([b(:,[1,2]); b(:,[2,3]); b(:,[3,4]); b(:,[4,1])],2),'rows');
L = sqrt(sum((p(e(:,1),:)-p(e(:,2),:)).^2,2));                             % Edge length

ne = size(e,1);
for k = 1:ne                                                               % Initial h - minimum neighbouring edge length
   Lk = L(k);
   if Lk<h(e(k,1)), h(e(k,1)) = Lk; end             
   if Lk<h(e(k,2)), h(e(k,2)) = Lk; end
end
h = min(h,hmax);

% Gradient limiting
tol = 1.0e-06;
while true                                                                 % Loop over the edges of the background mesh ensuring
   h_old = h;                                                              % that dh satisfies the dhmax tolerance
   for k = 1:ne                                                            % Loop over edges
      n1 = e(k,1);
      n2 = e(k,2);
      Lk = L(k);
      if h(n1)>h(n2)                                                       % Ensure grad(h)<=dhmax
         dh = (h(n1)-h(n2))/Lk;
         if dh>dhmax
            h(n1) = h(n2) + dhmax*Lk;
         end
      else
         dh = (h(n2)-h(n1))/Lk;
         if dh>dhmax
            h(n2) = h(n1) + dhmax*Lk;
         end
      end
   end
   if norm((h-h_old)./h,'inf')<tol                                         % Test convergence
      break
   end
end

% TRIANGULATE QUADTREE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if output
   fprintf('Triangulating quadtree\n');
end

if size(b,1)==1
   % Split box diagonally into 2 tri's
   t = [b(1),b(2),b(3); b(1),b(3),b(4)];
else

   % Get node-to-node connectivity
   % First column is column count per row
   % Max neighbours is 8 due to quadtree setup
   n2n = zeros(size(p,1),9);
   for k = 1:size(e,1)
      % Nodes in kth edge
      n1 = e(k,1);
      n2 = e(k,2);
      % Indexing
      n2n(n1,1) = n2n(n1,1)+1;                                             % Node 1
      n2n(n1,n2n(n1,1)+1) = n2;
      n2n(n2,1) = n2n(n2,1)+1;                                             % Node 2
      n2n(n2,n2n(n2,1)+1) = n1;
   end

   % Check for regular boxes with no mid-side nodes
   num = n2n(:,1)<=4;
   reg = all(num(b),2);

   % Split regular boxes diagonally into 2 tri's
   t = [b(reg,[1,2,3]); b(reg,[1,3,4])];
  
   if ~all(reg)
      
      % Use the N2N connectivity to directly triangulate the quadtree
      % nodes. Some additional nodes may be added at the centroids of some
      % boxes to facilitate triangulation. The triangluation is not
      % necessarily Delaunay, but should always be high quality and
      % symmetric where possible.

      b = b(~reg,:);                                                       % Boxes that still need to be dealt with
      nb = size(b,1);
      nt = size(t,1);
      
      % Alloc space
      t = [t; zeros(5*nb,3)];                                              % Has to be a least 5 times as many tri's as boxes
      nlist = zeros(512,1);                                                % Shouldn't ever be exceeded!

      for k = 1:nb

         % Corner nodes
         n1 = b(k,1); n2 = b(k,2);
         n3 = b(k,3); n4 = b(k,4);

         % Assemble node list for kth box in CCW order
         nlist(1) = n1;
         count = 1;
         next = 2;
         while true
            
            cn = nlist(next-1);

            % Find the closest node to CN travelling CCW around box
            old = inf;
            for j = 1:n2n(cn,1)
               nn = n2n(cn,j+1);
               dx = p(nn,1)-p(cn,1);
               dy = p(nn,2)-p(cn,2);
               if count==1                         % Edge 12
                  if (dx>0.0) && (dx<old)
                     old = dx;
                     tmp = nn;
                  end
               elseif count==2                     % Edge 23
                  if (dy>0.0) && (dy<old)
                     old = dy;
                     tmp = nn;
                  end
               elseif count==3                     % Edge 34
                  if (dx<0.0) && (abs(dx)<old)
                     old = abs(dx);
                     tmp = nn;
                  end
               else                                % Edge 41
                  if (dy<0.0) && (abs(dy)<old)
                     old = abs(dy);
                     tmp = nn;
                  end
               end

            end
            
            if tmp==n1                                                     % Back to start - Done!
               break
            elseif (count<4) && (tmp==b(k,count+1))                        % New edge
               count = count+1;
            end
            nlist(next) = tmp;
            next = next+1;
            
         end
         nnode = next-1;

         if (nt+nnode)>=size(t,1)                                          % Realloc memory on demand
            t = [t; zeros(nb,3)];
         end
         if (np+1)>=size(p,1)
            p = [p; zeros(nb,2)];
            h = [h; zeros(nb,1)];
         end
         
         % Triangulate box
         if nnode==4                                                       % Special treatment if no mid-side nodes
                                                                           % Split box diagonally into 2 tri's
            % New tri's
            t(nt+1,1) = n1;                     % t1
            t(nt+1,2) = n2;
            t(nt+1,3) = n3;
            t(nt+2,1) = n1;                     % t2
            t(nt+2,2) = n3;
            t(nt+2,3) = n4;
            
            % Update count
            nt = nt+2;
            
         elseif nnode==5                                                   % Special treatment if only 1 mid-side node
                                                                           % Split box into 3 tri's centred at mid-side node
            % Find the mid-side node
            j = 2;
            while j<=4
               if nlist(j)~=b(k,j)
                  break
               end
               j = j+1;
            end
           
            % Permute indexing so that the split is always between n1,n2
            if j==3
               n1 = b(k,2);   n2 = b(k,3);
               n3 = b(k,4);   n4 = b(k,1);
            elseif j==4
               n1 = b(k,3);   n2 = b(k,4);
               n3 = b(k,1);   n4 = b(k,2);
            elseif j==5
               n1 = b(k,4);   n2 = b(k,1);
               n3 = b(k,2);   n4 = b(k,3);
            end
            
            % New tri's
            t(nt+1,1) = n1;                     % t1
            t(nt+1,2) = nlist(j);
            t(nt+1,3) = n4;
            t(nt+2,1) = nlist(j);               % t2
            t(nt+2,2) = n2;
            t(nt+2,3) = n3;
            t(nt+3,1) = n4;                     % t3
            t(nt+3,2) = nlist(j);
            t(nt+3,3) = n3;

            % Update count
            nt = nt+3;
            
         else                                                              % Connect all mid-side nodes to an additional node
                                                                           % introduced at the centroid
            % New tri's
            xave = 0.0;
            yave = 0.0;
            have = 0.0;
            for j = 1:nnode-1
               jj = nlist(j);
               % New tri's
               t(nt+j,1) = jj;
               t(nt+j,2) = np+1;
               t(nt+j,3) = nlist(j+1);
               % Averaging
               xave = xave+p(jj,1);
               yave = yave+p(jj,2);
               have = have+h(jj);
            end
            jj = nlist(nnode);
            % Last tri
            t(nt+nnode,1) = jj;
            t(nt+nnode,2) = np+1;
            t(nt+nnode,3) = nlist(1);
            % New node
            p(np+1,1)   = (xave+p(jj,1)) /nnode;
            p(np+1,2)   = (yave+p(jj,2)) /nnode;
            h(np+1)     = (have+h(jj))   /nnode;

            % Update count
            nt = nt+nnode;
            np = np+1;

         end

      end
      p = p(1:np,:);
      h = h(1:np);
      t = t(1:nt,:);

   end

end

% Undo rotation
p = rotate(p,-theta);

end      % quadtree()


%% SUB-FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta = minrectangle(p)

% Find the rotation angle that must be applied to the 2D points in P so
% that the long axis of the minimum bounding rectangle is aligned with the 
% Y axis.

n = size(p,1);
if numel(p)~=2*n
   error('P must be an Nx2 array');
end

if n>2
   
   % Convex hull edge segments
   e = convhulln(p);
   
   % Keep convex points
   i = unique(e(:));
   p = p(i,:);
   
   % Re-index to keep E consistent
   j = zeros(size(p,1),1);
   j(i) = 1;
   j = cumsum(j);
   e = j(e);
   
   % Angles of hull segments
   dxy = p(e(:,2),:)-p(e(:,1),:);
   ang = atan2(dxy(:,2),dxy(:,1));            
   
   % Check all hull edge segments
   Aold = inf;
   for k = 1:size(e,1)
      % Rotate through -ang(k)
      pr = rotate(p,-ang(k));
      % Compute area of bounding rectangle and save if better
      dxy = max(pr,[],1)-min(pr,[],1);
      A = dxy(1)*dxy(2);
      if A<Aold
         Aold = A;
         theta = -ang(k);
      end
   end
   
   % Check result to ensure that the long axis is aligned with Y
   pr = rotate(p,theta);
   dxy = max(pr,[],1)-min(pr,[],1);
   if dxy(1)>dxy(2)
      % Need to flip XY
      theta = theta+0.5*pi;
   end
   
else
   % 1 or 2 points, degenerate bounding rectangle in either case
   theta = 0.0;
end

end      % minrectangle()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = rotate(p,theta)

% Rotate the 2D points in P through the angle THETA (radians).

stheta = sin(theta);
ctheta = cos(theta);

p = p*[ctheta, stheta; -stheta, ctheta];

end      % rotate()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = userhfun(x,y,fun,args,hmax,xymin,xymax)

% Evaluate user defined size function.

if ~isempty(fun)
   h = feval(fun,x,y,args{:});
   if size(h)~=size(x)
      error('Incorrect user defined size function. SIZE(H) must equal SIZE(X).');
   end
else
   h = inf*ones(size(x));
end
h = min(h,hmax);

% Limit to domain
out = (x>xymax(1))|(x<xymin(1))|(y>xymax(2))|(y<xymin(2));
h(out) = inf;

if any(h<=0.0)
   error('Incorrect user defined size function. H must be positive.');
end

end      % userhfun()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hmax,edgeh,fun,args] = gethdata(hdata)

% Check the user defined size functions

d_hmax  = inf;
d_edgeh = [];
d_fun   = '';
d_args  = {};

if ~isempty(hdata)
   if ~isstruct(hdata)
      error('HDATA must be a structure');
   end
   if numel(hdata)~=1
      error('HDATA cannot be an array of structures');
   end
   fields = fieldnames(hdata);
   names = {'hmax','edgeh','fun','args'};
   for k = 1:length(fields)
      if ~any(strcmp(fields{k},names))
         error('Invalid field in HDATA');
      end
   end
   if isfield(hdata,'hmax')
      if (numel(hdata.hmax)~=1) || (hdata.hmax<=0)
         error('HDATA.HMAX must be a positive scalar');
      else
         hmax = hdata.hmax;
      end
   else
      hmax = d_hmax;
   end
   if isfield(hdata,'edgeh')
      if (numel(hdata.edgeh)~=2*size(hdata.edgeh,1)) || any(hdata.edgeh(:)<0)
         error('HDATA.EDGEH must be a positive Kx2 array');
      else
         edgeh = hdata.edgeh;
      end
   else
      edgeh = d_edgeh;
   end
   if isfield(hdata,'fun')
      if ~ischar(hdata.fun) && ~isa(hdata.fun,'function_handle')
         error('HDATA.FUN must be a function name or a function handle');
      else
         fun = hdata.fun;
      end
   else
      fun = d_fun;
   end
   if isfield(hdata,'args')
      if ~iscell(hdata.args)
         error(['HDATA.ARGS must be a cell array of additional' ...
            'inputs for HDATA.FUN']);
      else
         args = hdata.args;
      end
   else
      args = d_args;
   end
else
   hmax  = d_hmax;
   edgeh = d_edgeh;
   fun   = d_fun;
   args  = d_args;
end

end      % gethdata()
