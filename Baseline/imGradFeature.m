% imgradfeature transforms am image to a feature representation which consists of intensity, forward horizontal diference, forward vertical diference, backward horizontal diference, and backward vertical diference
%
% Y = imGradFeature(X)
%Output parameters:
% Y(:,:,:,1): intensity
% Y(:,:,:,2): forward horizontal difference
% Y(:,:,:,3): forward vertical difference
% Y(:,:,:,4): backward horizontal difference
% Y(:,:,:,5): backward vertical difference
%
%
%Input parameters:
% X: input image
%
%
%Example:
% X = double(imread('img.png'));
% F = imGradFeature(X);
%
%
%Version: 20121212

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified Poisson                                         %
%                                                          %
% Copyright (C) 2012 Masayuki Tanaka. All rights reserved. %
%                    mtanaka@ctrl.titech.ac.jp             %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = imGradFeature(X)

s1 = size(X,1);
s2 = size(X,2);
s3 = 1;
if( ndims(X) >= 3 )
 s3 = size(X,3);
end

Y = zeros(s1,s2,s3,5);

Kh = [ 0,-1, 1 ];
Kv = [ 0;-1; 1 ];

Y(:,:,:,1) = X;
Y(:,:,:,2) = imfilter(X,Kh,'replicate');
Y(:,:,:,3) = imfilter(X,Kv,'replicate');
Y(:,:,:,4) = circshift(Y(:,:,:,2),[0,1]);
Y(:,:,:,5) = circshift(Y(:,:,:,3),[1,0]);

