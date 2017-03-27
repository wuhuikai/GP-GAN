function [p,t,fnum,stats] = meshfaces(node,edge,face,hdata,options)

%  MESHFACES: 2D unstructured mesh generation for polygonal geometry.
%
% A 2D unstructured triangular mesh is generated based on a piecewise-
% linear geometry input. An arbitrary number of polygonal faces can be 
% specified, and each face can contain an arbitrary number of cavities. An 
% iterative method is implemented to optimise mesh quality. 
%
% If you wish to mesh a single face, use MESH2D instead!
%
%  [p,t,fnum] = meshfaces(node,edge,face,hdata,options);
%
% OUTPUTS
%
%  P     = Nx2 array of nodal XY co-ordinates.
%  T     = Mx3 array of triangles as indicies into P, defined with a
%          counter-clockwise node ordering.
%  FNUM  = Mx1 array of face numbers for each triangle in T.
%
% INPUTS
%  
% Blank arguments can be passed using the empty placeholder "[]".
%
% NODE defines the XY co-ordinates of the geometry vertices. The vertices
% for all faces must be specified:
%
%  NODE = [x1,y1; x2,y2; etc], xy geometry vertices specified in any order.
%
% EDGE defines the connectivity between the points in NODE as a list of
% edges. Edges for all faces must be specified:
%
%  EDGE = [n1 n2; n2 n3; etc], connectivity between nodes to form
%                              geometry edges.
%
% FACE defines the edges included in each geometry face. Each face is a 
% vector of edge numbers, stored in a cell array:
%
%  FACE{1} = [e11,e12,etc]
%  FACE{2} = [e21,e22,etc]
%
% HDATA is a structure containing user defined element size information. 
% HDATA can include the following fields:
%
%  hdata.hmax  = h0;                   Max allowable global element size.
%  hdata.edgeh = [e1,h1; e2,h2; etc];  Element size on specified geometry 
%                                      edges.
%  hdata.fun   = 'fun' or @fun;        User defined size function.
%  hdata.args  = {arg1, arg2, etc};    Additional arguments for HDATA.FUN.
%
% Calls to user specified functions must accept vectorised input of the 
% form H = FUN(X,Y,ARGS{:}), where X,Y are the xy coordinates where the
% element size will be evaluated and ARGS are optional additional arguments 
% as passed by HDATA.ARGS.
%
% An automatic size function is always generated to ensure that the
% geometry is adequately resolved. The overall size function is the minimum
% of the user specified and automatic functions.
%
% OPTIONS is a structure array that allows some of the "tuning" parameters
% used in the solver to be modified:
%
%   options.mlim   : The convergence tolerance. The maximum percentage 
%                    change in edge length per iteration must be less than 
%                    MLIM { 0.02, 2.0% }. 
%   options.maxit  : The maximum allowable number of iterations { 20 }.
%   options.dhmax  : The maximum allowable (relative) gradient in the size 
%                    function { 0.3, 30.0% }.
%   options.output : Displays the mesh and the mesh statistics upon
%                    completion { TRUE }.
%
% EXAMPLE:
%
%   meshdemo                  % Will run the standard demos
%   mesh_collection(n)        % Will run some additional demos
%
% See also MESH2D, REFINE, SMOOTHMESH, DELAUNAYN

% STATS is an undocumented output used in debugging. Returns the algorithm 
% statistics usually printed to screen as a structure.

%   Darren Engwirda : 2005-09
%   Email           : d_engwirda@hotmail.com
%   Last updated    : 10/10/2009 with MATLAB 7.0 (Mesh2d v2.4)
%
% Please email me any un-meshable geometries, meshing benchmarks or
% suggestions!

ts = cputime;

% Catch any initalisation errors
try                                                   

   % I/O error checks
   if (nargin<5)
      options = [];
      if (nargin<4)
         hdata = [];
         if (nargin<3)
            face = [];
            if (nargin<2)
               edge = [];
               if (nargin<1)
                  error('Wrong number of inputs');
               end
            end
         end
      end
   elseif (nargin>5)
      error('Wrong number of inputs');
   end
   if (nargout>4)
      error('Wrong number of outputs');
   end

   % Get user options
   options = getoptions(options);
   
   % Check geometry and attempt to repair bad geometry
   if options.output
      fprintf('Checking Geometry\n');
   end
   [node,edge,face,hdata] = checkgeometry(node,edge,face,hdata);      
   
catch
   % Close waitbar on error
   close(wbar);
   rethrow(lasterror);
end

% Quadtree decomposition
%  PH    : Background mesh nodes
%  TH    : Background mesh triangles
%  HH    : Size function value at PH
tic
[qtree.p,qtree.t,qtree.h] = quadtree(node,edge,hdata,options.dhmax,options.output);
t_quad = toc;

% Discretise edges
pbnd = boundarynodes(qtree.p,qtree.t,qtree.h,node,edge,options.output);

% Mesh each face separately
p = []; t = []; fnum = [];
for k = 1:length(face)
   
   % Mesh kth polygon
   [pnew,tnew] = meshpoly(node,edge(face{k},:),qtree,pbnd,options);

   % Add to global lists
   t = [t; tnew+size(p,1)];
   p = [p; pnew];
   fnum = [fnum; k*ones(size(tnew,1),1)];
   
end

% Ensure consistent, CCW ordered triangulation
[p,t,fnum,fnum] = fixmesh(p,t,[],fnum);   

% Element quality
q = quality(p,t);

% Method statistics
stats = struct('Time',cputime-ts,'Triangles',size(t,1), ...
                  'Nodes',size(p,1),'Mean_quality',mean(q),'Min_quality',min(q));
               
if options.output
   %figure('Name','Mesh')
   %plot(p(:,1),p(:,2),'b.','markersize',1)
   %hold on;
   % Colour mesh for each face
%    col = ['b','r','g','o','m'];
%    for k = 1:length(face)
%       colk = mod(k,length(col));
%       if (colk==0)
%          colk = length(col);
%       end
%       patch('faces',t(fnum==k,:),'vertices',p,'facecolor','w','edgecolor',col(colk));
%    end
%    patch('faces',edge,'vertices',node,'facecolor','none','edgecolor','k')
%    % Highlisght low q triangles in debug mode
%    if options.debug
%       pc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3.0;
%       plot(pc(q<0.5,1),pc(q<0.5,2),'r.')
%    end
%    axis equal off;
%    disp(stats);
end

end      % meshfaces()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = boundarynodes(ph,th,hh,node,edge,output)

% Discretise the geometry based on the edge size requirements interpolated 
% from the background mesh.

p = node;
e = edge;
i = tsearchn([ph(:,1),ph(:,2)],th,[p(:,1),p(:,2)]);               
h = tinterp(ph,th,hh,p,i);

if output
   fprintf('Placing Boundary Nodes\n');
end
iter = 1;
while true
   
   % Edge length
   dxy = p(e(:,2),:)-p(e(:,1),:);
   L = sqrt(sum(dxy.^2,2));
   % Size function on edges                                           
   he = 0.5*(h(e(:,1))+h(e(:,2)));
   % Split long edges
   ratio = L./he;
   split = (ratio>=1.5);
   if any(split)
      % Split edge at midpoint
      n1 = e(split,1);
      n2 = e(split,2);
      pm = 0.5*(p(n1,:)+p(n2,:));
      n3 = (1:size(pm,1))' + size(p,1);
      % New lists
      e(split,:) = [n1,n3];
      e = [e; n3,n2];
      p = [p; pm];
      % Size function at new nodes
      i = mytsearch(ph(:,1),ph(:,2),th,pm(:,1),pm(:,2));               
      h = [h; tinterp(ph,th,hh,pm,i)];
   else
      break
   end
   iter = iter+1;
end

% Node-to-edge connectivity matrix
ne = size(e,1);
S = sparse(e(:),[1:ne,1:ne],[-ones(ne,1); ones(ne,1)],size(p,1),ne);

% Smooth bounday nodes
if output
   fprintf('Smoothing Boundaries\n');
end
del = 0.0;
tol = 0.02;
maxit = 50;
i = zeros(size(p,1),1);
for iter = 1:maxit 

   delold = del;
   
   % Spring based smoothing
   F = he./L-1.0;
   F = S*(dxy.*[F,F]);
   F(1:size(node,1),:) = 0.0;   
   p = p+0.2*F;

   % Convergence
   dxy = p(e(:,2),:)-p(e(:,1),:);
   Lnew = sqrt(sum(dxy.^2,2));
   del = norm((Lnew-L)./Lnew,'inf');
   if (del<tol)
      break;
   else
      if (iter==maxit)
         disp('WARNING: Boundary smoothing did not converge.');
      end
   end
   L = Lnew;
   
   if (del>delold)
      % Interpolate required size at new P
      i = mytsearch(ph(:,1),ph(:,2),th,p(:,1),p(:,2),i);
      h = tinterp(ph,th,hh,p,i);
      he = 0.5*(h(e(:,1))+h(e(:,2)));
   end
   
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function options = getoptions(options)

% Extract the user defined options

% Defaults
d_mlim   = 0.02;
d_maxit  = 20;
d_dhmax  = 0.3;
d_output = true;

if ~isempty(options)
   if ~isstruct(options)
      error('OPTIONS must be a structure array');
   end
   if numel(options)~=1
      error('Options cannot be an array of structures');
   end
   fields = fieldnames(options);
   names = {'mlim','maxit','dhmax','output'};
   for k = 1:length(fields)
      if strcmp(fields{k},names)
         error('Invalid field in OPTIONS');
      end
   end
   if isfield(options,'mlim')                                              % Movement tolerance
      checkposscalar(options.mlim,'options.mlim');
   else
      options.mlim = d_mlim;
   end
   if isfield(options,'maxit')                                             % Maximum iterations
      options.maxit = round(checkposscalar(options.maxit,'options.maxit'));
   else
      options.maxit = d_maxit;
   end
   if isfield(options,'dhmax')                                             % Size function gradient limit
      checkposscalar(options.dhmax,'options.dhmax');
   else
      options.dhmax = d_dhmax;
   end
   if isfield(options,'output')                                            % Output on/off
      checklogicalscalar(options.output,'options.output');
   else
      options.output = d_output;
   end
else                                                                       % Default values
   options.mlim   = d_mlim;
   options.maxit  = d_maxit;
   options.dhmax  = d_dhmax;
   options.output = d_output;
end
options.debug = false;
options.output = false;

end      % getoptions()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var = checkposscalar(var,name)

% Helper function to check if var is a positive scalar.

if var<0 || any(size(var)>1)
   error([name,' must be a positive scalar']);
end

end      % checkposscalar()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var = checklogicalscalar(var,name)

% Helper function to check if var is a logical scalar.

if ~islogical(var) || any(size(var)>1)
   error([name,' must be a logical scalar']);
end

end      % checklogicalscalar()
