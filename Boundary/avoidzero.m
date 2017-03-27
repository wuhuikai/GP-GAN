% avoidzero to avoid zero
%
% Y = avoidzero( X, smallpositive )
%Output parameters:
% Y: zero-avoided data
%
%
%Input parameters:
% X: the input data
% smallpositive : small positive value to avoid divided-by-zero
%
%
%Version: 20160622

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified Poisson                                         %
%                                                          %
% Copyright (C) 2016 Masayuki Tanaka. All rights reserved. %
%                    mtanaka@ctrl.titech.ac.jp             %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = avoidzero( X, smallpositive )

Y = X;
Y(find( X >= 0 & X < smallpositive )) = smallpositive;
Y(find( X < 0 & X > -smallpositive )) = -smallpositive;

end
