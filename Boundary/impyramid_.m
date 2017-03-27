function B = impyramid(A, direction)
%IMPYRAMID Image pyramid reduction and expansion
%    B = IMPYRAMID(A, DIRECTION) computes a Gaussian pyramid reduction or
%    expansion of A by one level.  DIRECTION can be 'reduce' or 'expand'.  
% 
%    If A is M-by-N and DIRECTION is 'reduce', then the size of B is
%    ceil(M/2)-by-ceil(N/2).  If DIRECTION is 'expand', then the size of B is
%    (2*M-1)-by-(2*N-1). 
%
%    Reduction and expansion take place only in the first two dimensions.  For
%    example, if A is 100-by-100-by-3 and DIRECTION is 'reduce', then B is
%    50-by-50-by-3.
% 
%    Class support
%    -------------
%    A can be any numeric class except uint64 or int64, or it can be 
%    logical.  The class of B is the same as the class of A.
% 
%    Example
%    -------
%    Compute a four-level multiresolution pyramid of the cameraman image.
% 
%       I0 = imread('cameraman.tif');
%       I1 = impyramid(I0, 'reduce');
%       I2 = impyramid(I1, 'reduce');
%       I3 = impyramid(I2, 'reduce');
% 
%       imshow(I0)
%       figure, imshow(I1)
%       figure, imshow(I2)
%       figure, imshow(I3)
%
%   See also IMRESIZE.

%   Copyright 2011 The MathWorks, Inc.

% References:
% Burt and Adelson, "The Laplacian Pyramid as a Compact Image Code," IEEE
% Transactions on Communications, vol. COM-31, no. 4, April 1983,
% pp. 532-540.
%
% Burt, "Fast Filter Transforms for Image Processing," Computer Graphics and
% Image Processing, vol. 16, 1981, pp. 20-51.

validateattributes(A, {'numeric', 'logical'}, {}, mfilename, 'A', 1);
direction = validatestring(direction, {'reduce', 'expand'}, ...
    mfilename, 'DIRECTION', 2);

M = size(A,1);
N = size(A,2);

if strcmp(direction, 'reduce')
    scaleFactor = 0.5;
    outputSize = ceil([M N]/2);
    kernel = makePiecewiseConstantFunction( ...
        [3.5 2.5      1.5    0.5    -0.5   -1.5    -Inf], ...
        [0.0 0.0625   0.25   0.375   0.25   0.0625  0.0]);
    kernelWidth = 5;
    
else
    scaleFactor = 2;
    outputSize = 2*[M N];
    kernel = makePiecewiseConstantFunction( ...
        [1.25   0.75    0.25   -0.25   -0.75   -1.25   -Inf], ...
        [0.0    0.125   0.5     0.75    0.5    0.125    0.0]);
    kernelWidth = 3;
end

B = imresize(A, scaleFactor, {kernel, kernelWidth}, ...
    'OutputSize', outputSize, 'Antialiasing', false);

end

function fun = makePiecewiseConstantFunction(breakPoints, values)
% Constructs a piecewise constant function and returns a handle to it.
%
% breakPoints and values have to be vectors with the same number of
% elements.
%
% The elements in breakPoints have to be monotonically decreasing.
% 
% fun(x) returns values(1) if x >= breakPoints(1)
%
% else fun(x) returns values(2) if x >= breakPoints(2)
%
% else fun(x) returns values(3) if x >= breakPoints(3)
%
% etc.
%
% If x is less than breakPoint(end), then fun returns 0.
%
% If x is an array, then fun operates elementwise on x and returns an array
% of the same size.

%iptassert(all(diff(breakPoints) < 0), ...
%          'images:impyramid:badBreakPointList')

fun = @piecewiseConstantFunction;

    function y = piecewiseConstantFunction(x)
        y = zeros(size(x));
        for k = 1:numel(x)
            yy = 0;
            xx = x(k);
            for p = 1:numel(breakPoints)
                if xx >= breakPoints(p)
                    yy = values(p);
                    break;
                end
            end
            y(k) = yy;
        end
    end
end
