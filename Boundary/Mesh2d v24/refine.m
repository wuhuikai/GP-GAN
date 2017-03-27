function [p,t,f] = refine(p,t,ti,f)

%  REFINE: Refine triangular meshes.
%
% Quadtree triangle refinement is performed, with each triangle split into
% four sub-triangles. The new triangles are created by joining nodes
% introduced at the edge midpoints. The refinement is "quality" preserving,
% with the aspect ratio of the sub-triangles being equal to that of the
% parent.
%
%  UNIFORM REFINEMENT:
%
%  [p,t] = refine(p,t);
%
%  p : Nx2 array of nodal XY coordinates, [x1,y1; x2,y2; etc]
%  t : Mx3 array of triangles as indices, [n11,n12,n13; n21,n22,n23; etc]
%
%  NON-UNIFORM REFINEMENT:
%
% Non-uniform refinement can also be performed by specifying which
% triangles are to be refined. Quadtree refinement is performed on
% specified triangles. Neighbouring triangles are also refined in order to
% preserve mesh compatibility. These triangles are refined using
% bi-section.
%
%  [p,t] = refine(p,t,ti);
%
%  ti : Mx1 logical array, with Ti(k) = TRUE if kth triangle is to be
%       refined
%
% Functions defined on the nodes in P can also be refined using linear
% interpolation through an extra input:
%
%  [p,t,f] = refine(p,t,ti,f);
%
%  f : NxK array of nodal function values. Each column in F corresponds to
%      a dependent function, F(:,1) = F1(P), F(:,2) = F2(P), etc.
%
% It is often useful to smooth the refined mesh using SMOOTHMESH. Generally
% this will improve element quality.
%
% Example:
%
%   [p,t] = refine(p,t,ti);
%
% See also SMOOTHMESH, MESH2D

%   Darren Engwirda : 2007
%   Email           : d_engwirda@hotmail.com
%   Last updated    : 26/07/2007 with MATLAB 7.0

if nargin<=4
   gotF = false;
   if (nargin<=2) || isempty(ti)                                           % Uniform refinement
      ti = true(size(t,1),1);
      if nargin<2
         error('Wrong number of inputs');
      end
   end
else
   gotF = true;
   if nargin>4
      error('Wrong number of inputs');
   end
end
if (gotF&&(nargout>3)) || (~gotF&&(nargout>2))
   error('Wrong number of outputs');
end
if numel(ti)~=size(t,1)
   error('Ti must be an Mx1 array');
end
if gotF && ((size(f,1)~=size(p,1)) || (ndims(f)>2))
   error('F must be an NxK array');
end

% Ensure we start with a consistent mesh
if gotF
   [p,t,junk,f] = fixmesh(p,t,[],f);
else
   [p,t] = fixmesh(p,t);
end

% Edge connectivity
numt = size(t,1);
vect = 1:numt;
e = [t(:,[1,2]); t(:,[2,3]); t(:,[3,1])];                                  % Edges - not unique
[e,j,j] = unique(sort(e,2),'rows');                                        % Unique edges
te = [j(vect), j(vect+numt), j(vect+2*numt)];                              % Unique edges in each triangle

split = false(size(e,1),1);
split(te(ti,:)) = true;                                                    % True for edges to be split

% Flag tri's to be split
nsplit = length(find(split));
while true
   split3 = sum(double(split(te)),2)>=2;                                   % True for tri's where we will split 3 edges
   split(te(split3,:)) = true;                                             % Update split - the split2 case was turned into
   new = length(find(split))-nsplit;                                       % a split3 case, setting new edges to be split
   if new==0
      break
   end
   nsplit = nsplit+new;
end
split1 = sum(double(split(te)),2)==1;                                      % True for tri's where we will split 1 edge

% New nodes
np = size(p,1);
pm = 0.5*(p(e(split,1),:)+p(e(split,2),:));                                % Split edge midpoints
p = [p; pm];

% Map E(SPLIT) to index PM
i = zeros(size(e,1),1);
i(split) = (1:nsplit)'+np;

% New tri's in the split3 case
tnew = t(~(split1|split3),:);
if any(split3)
   n1 = t(split3,1);
   n2 = t(split3,2);
   n3 = t(split3,3);
   n4 = i(te(split3,1));
   n5 = i(te(split3,2));
   n6 = i(te(split3,3));
   tnew = [tnew; n1,n4,n6; n4,n2,n5; n5,n3,n6; n4,n5,n6];
end

% New tri's in the split1 case
if any(split1)
   [row,col] = find(split(te(split1,:)));                                  % Find split edges in tri's

   N1 = col;                                                               % Transform so that the split is always between n1 & n2
   N2 = col+1;
   N3 = col+2;
   N2(N2>3) = N2(N2>3)-3;
   N3(N3>3) = N3(N3>3)-3;

   n1 = 0*N1;
   n2 = n1;
   n3 = n1;
   n4 = n1;

   split1 = find(split1);
   split1 = split1(row);
   for k = 1:length(col)
      n1(k) = t(split1(k),N1(k));
      n2(k) = t(split1(k),N2(k));
      n3(k) = t(split1(k),N3(k));
      n4(k) = i(te(split1(k),col(k)));
   end
   tnew = [tnew; n1,n4,n3; n4,n2,n3];
end
t = tnew;

% Linear interpolation to new nodes
if gotF
   f = [f; 0.5*(f(e(split,1),:)+f(e(split,2),:))];
end

end      % refine()
