% buildModPoissonParam provides a parameters which is used in ModPoisson.
%
% param = buildModPoissonParam( s )
%Output parameters:
% param: the parameters to be used in ModPoisson
%
%
%Input parameters:
% s: [s1, s2] size of output image
%
%
%Example:
% param = buildModPoissonParam( [256, 256] );
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
function param = buildModPoissonParam( s1, s2 )

if( nargin == 1 )
 s = s1;
else
 s = [s1 s2];
end

K=zeros(s(1)*2,s(2)*2);
K(1,1)=4;
K(1,2)=-1;
K(2,1)=-1;
K(s(1)*2,1)=-1;
K(1,s(2)*2)=-1;
param = fft2(K);
param = real(param(1:s(1),1:s(2)));
