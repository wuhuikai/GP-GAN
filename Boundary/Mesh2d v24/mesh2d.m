function [p,t,stats] = mesh2d(node,edge,hdata,options)

%  MESH2D: 2D unstructured mesh generation for a polygon.
%
% A 2D unstructured triangular mesh is generated based on a piecewise-
% linear geometry input. The polygon can contain an arbitrary number of 
% cavities. An iterative method is implemented to optimise mesh quality. 
%
% If you wish to mesh multiple connected faces, use MESHFACES instead!
%
% OUTPUTS
%
%  P     = Nx2 array of nodal XY co-ordinates.
%  T     = Mx3 array of triangles as indicies into P, defined with a
%          counter-clockwise node ordering.
%
% SHORT SYNTAX:
%
%  [p,t] = mesh2d(node);
%
% NODE defines the geometry nodes as an Nx2 array:
%
%  node  = [x1 y1; x2 y2; etc], geometry nodes specified in consectutive
%                               order, such that NODE(2,:) is joined with
%                               NODE(1,:) etc.
%
% An element size function is automatically generated based on the 
% complexity of the geometry. Generally this produces meshes with the 
% fewest number of triangles.
%
% LONG SYNTAX:
%
%  [p,t] = mesh2d(node,edge,hdata,options);
%
% Blank arguments can be passed using the empty placeholder "[]".
%
% EDGE defines the connectivity between the points in NODE as a list of
% edges:
%
%   edge = [n1 n2; n2 n3; etc], connectivity between nodes to form
%                               geometry edges. If EDGE is specified it is
%                               not required that NODE be consectutive.
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
% See also MESHFACES, REFINE, SMOOTHMESH, DELAUNAYN

% STATS is an undocumented output used in debugging. Returns the algorithm 
% statistics usually printed to screen as a structure.

%   Darren Engwirda : 2005-09
%   Email           : d_engwirda@hotmail.com
%   Last updated    : 10/10/2009 with MATLAB 7.0 (Mesh2d v2.4)
%
% Please email me any un-meshable geometries, meshing benchmarks or
% suggestions!

if (nargin<4)
   options = [];
   if (nargin<3)
      hdata = [];
      if (nargin<2)
         edge = [];
      end
   end
end

% Assume 1 face containing all edges
[p,t,junk,stats] = meshfaces(node,edge,[],hdata,options);

end      % mesh2d()
