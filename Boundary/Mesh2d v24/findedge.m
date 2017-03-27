function enum = findedge(p,node,edge,TOL)

%  FINDEDGE: Locate points on edges.
%
% Determine which edges a series of points lie on in a 2D plane. 
%
%  i = findedge(p,node,edge,tol);
%
% INPUTS
%
%  P     : An Nx2 array of xy co-ordinates of points to be checked.
%  NODE  : An Kx2 array of xy co-ordinates of edge endpoints.
%  EDGE  : An Mx2 array of edges, specified as connections between the 
%          vertices in NODE: [n1 n2; n3 n4; etc]. 
%  TOL   : Tolerance used when testing points.
%
% OUTPUTS
%
%  I     : Nx1 array of edge numbers, corresponding to the edge that each
%          node lies on. Nodes that do not lie on any edges are assigned 0.
%
% See also INPOLYGON

%   Darren Engwirda: 2005-2007
%   Email          : d_engwirda@hotmail.com
%   Last updated   : 25/11/2007 with MATLAB 7.0
%
% Problems or suggestions? Email me.

%% ERROR CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
   TOL = 1.0e-12;
   if nargin<3
      edge = [];
      if nargin<2
         error('Insufficient inputs');
      end
   end
end
nnode = size(node,1);
if isempty(edge)                                                           % Build edge if not passed
   edge = [(1:nnode-1)' (2:nnode)'; nnode 1];
end
if size(p,2)~=2
   error('P must be an Nx2 array.');
end
if size(node,2)~=2
   error('NODE must be an Mx2 array.');
end
if size(edge,2)~=2
   error('EDGE must be an Mx2 array.');
end
if max(edge(:))>nnode || any(edge(:)<1)
   error('Invalid EDGE.');
end

%% PRE-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n  = size(p,1);
nc = size(edge,1);

% Choose the direction with the biggest range as the "y-coordinate" for the
% test. This should ensure that the sorting is done along the best
% direction for long and skinny problems wrt either the x or y axes.
dxy = max(p,[],1)-min(p,[],1);
if dxy(1)>dxy(2)
   % Flip co-ords if x range is bigger
   p = p(:,[2,1]);
   node = node(:,[2,1]);
end
tol = TOL*min(dxy);

% Sort test points by y-value
[y,i] = sort(p(:,2));
x = p(i,1);

%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enum = zeros(size(p,1),1);
for k = 1:nc         % Loop through edges

   % Nodes in current edge
   n1 = edge(k,1);
   n2 = edge(k,2);

   % Endpoints - sorted so that [x1,y1] & [x2,y2] has y1<=y2
   y1 = node(n1,2);
   y2 = node(n2,2);
   if y1<y2
      x1 = node(n1,1);
      x2 = node(n2,1);
   else
      yt = y1;
      y1 = y2;
      y2 = yt;
      x1 = node(n2,1);
      x2 = node(n1,1);
   end

   % Binary search to find first point with y<=y1 for current edge
   if y(1)>=y1
      start = 1;
   elseif y(n)<y1
      start = n+1;       
   else
      lower = 1;
      upper = n;
      for j = 1:n
         start = round(0.5*(lower+upper));
         if y(start)<y1
            lower = start;
         elseif y(start-1)<y1
            break;
         else
            upper = start;
         end
      end
   end

   % Loop through points
   for j = start:n
      % Check the bounding-box for the edge before doing the intersection
      % test. Take shortcuts wherever possible!
      Y = y(j);   % Do the array look-up once & make a temp scalar
      if Y<=y2
         
         % Check if we're "on" the edge
         X = x(j);
         if  (abs((y2-Y)*(x1-X)-(y1-Y)*(x2-X))<tol);
            enum(j) = k;
         end
         
      else
         % Due to the sorting, no points with >y
         % value need to be checked
         break
      end
   end

end

% Re-index to undo the sorting
enum(i) = enum;

end      % findedge()
