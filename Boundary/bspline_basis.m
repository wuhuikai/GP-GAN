function [value]=bspline_basis(x)
    
  value=(x+1).*(x>=-1 & x<0)+(1-x).*(x>=0 & x<=1);


end