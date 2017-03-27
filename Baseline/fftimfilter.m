
%the file was created by 
function res = fftimfilter(a,b)
%FFTIMFILTER imfilter via fft
%   Assumes size(a)==size(b)

sm = size(a);
res = conv2olam(a,b);
res = res(floor(sm(1)/2):end-ceil(sm(1)/2), floor(sm(2)/2):end-ceil(sm(2)/2));

end

