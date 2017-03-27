% the code is created according to Farbman's code
function [res]=PoissonBlending(trg,src,mask)
    h = fspecial('laplacian', 0);
    chi = imfilter(double(mask),h);
    chi(chi<0) = 0;
    chi(chi>0) = 1;
    erf = trg - src;
    res = zeros(size(erf));
    for i=1:3    
        sr = src(:,:,i);
        tr = trg(:,:,i);       
        a = erf(:,:,i);
        a(~chi) =  0;
        Ierf = LaplacianDirichlet(a,mask);
        temp = Ierf + sr;  
        tr(mask) = temp(mask);
        res(:,:,i) = tr;
    end   
end