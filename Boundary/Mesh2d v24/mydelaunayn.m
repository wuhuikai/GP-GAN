function t = mydelaunayn(p)

% My version of the MATLAB delaunayn function that attempts to deal with
% some of the compatibility and robustness problems.
%
% Darren Engwirda - 2007

% Translate to the origin and scale the min xy range onto [-1,1]
% This is absolutely critical to avoid precision issues for large problems!
maxxy = max(p);
minxy = min(p);
p(:,1) = p(:,1)-0.5*(minxy(1)+maxxy(1));
p(:,2) = p(:,2)-0.5*(minxy(2)+maxxy(2));
p = p/(0.5*min(maxxy-minxy));

try
   % Use the default settings. This will be 'Joggled Input' prior to
   % MATLAB 7.0 or 'Triangulated Output' in newer versions
   t = delaunayn(p);
catch
   if ~(str2double(version('-release'))<=13)
      % 'Triangulated Output' is generally more accurate, but less robust.
      % If Qhull crashes with 'Triangulated Output', redo with 'Joggled
      % Input'
      t = delaunayn(p,{'QJ','Pp'});
   end
end

end      % mydelaunayn()

% % !! If you copy "qhullmx.dll" into the folder containing "MyDelaunayn.m"
% % uncomment the following code to get a much faster implementation of the
% % call to Qhull. This is based on the standard MATLAB "delaunayn.m"
% % function, but optimised for 2D input.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function t = delaunayn(x,options)
% 
% % 2D optimised Qhull call.
% 
% if nargin < 1
%    error('Needs at least 1 input.');
% end
% if isempty(x)
%    t = []; 
%    return 
% end
% 
% [m,n] = size(x);
% if m < n+1,
%    error('Not enough unique points to do tessellation.');
% end
% if any(isinf(x(:)) | isnan(x(:)))
%    error('Data containing Inf or NaN cannot be tessellated.');
% end
% if m == n+1
%    t = 1:n+1;
%    return;
% end
% 
% % Deal with options
% if n >= 4
%    opt = 'Qt Qbb Qc Qx';
% else
%    opt = 'Qt Qbb Qc';
% end
% if ( nargin > 1 && ~isempty(options) )
%    if ~iscellstr(options)
%       error('OPTIONS should be cell array of strings.');
%    end
%    sp = {' '};
%    c = strcat(options,sp);
%    opt = cat(2,c{:});
% end
% 
% % Call Qhull to do the work
% t = qhullmx(x', 'd ', opt);
% 
% % Try to get rid of zero volume simplices. They are generated
% % because of the fuzzy jiggling.
% seps = eps^(4/5)*max(abs(x(:)));
% 
% % Triangle area
% d12 = x(t(:,2),:)-x(t(:,1),:);
% d13 = x(t(:,3),:)-x(t(:,1),:);
% A = (d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1));
% 
% t = t(abs(A)>seps,:);
% 
% end       % delaunauyn()
