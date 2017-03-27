function ahat = evalf( a, h1, h2, g )
%EVALF Evaluate convolution pyramid
%   Evaluate convolution pyramid on input a using 
%   filter set specified by: h1,h2,g.
%   For more details: 
%   http://www.cs.huji.ac.il/labs/cglab/projects/convpyr

[h,w] = size(a);
maxLevel = ceil(log2(max(h,w)));
fs = size(h1,1);

% Forward transform (analysis)
pyr{1} = padarray(a, [fs fs]);
for i=2:maxLevel
    
    down = imfilter(pyr{i-1},h1,0);
    down = down(1:2:end,1:2:end);
    
    down = padarray(down, [fs fs]);
    pyr{i} = down;
    
end

% Backward transform (synthesis)
fpyr{maxLevel} = imfilter(pyr{maxLevel},g,0);
for i=maxLevel-1:-1:1
    
    rd = fpyr{i+1};
    rd = rd(1+fs:end-fs, 1+fs:end-fs);
    
    up = zeros(size(pyr{i}));
    up(1:2:end,1:2:end) = rd;

    fpyr{i} = imfilter(up,h2,0) + imfilter(pyr{i},g,0);
end

ahat = fpyr{1};
ahat = ahat(1+fs:end-fs, 1+fs:end-fs);

end

