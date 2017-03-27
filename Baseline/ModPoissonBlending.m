%the code is created according to Masayuki Tanaka's code
function [res]=ModPoissonBlending(ftrg,src,mask)
    targetFeature = imGradFeature(ftrg);
    sourceFeature = imGradFeature(src);
    mask = repmat(double(mask), [1, 1, 3, 5]);
    targetFeature = targetFeature .* (1-mask) + sourceFeature .* mask;
    param = buildModPoissonParam(size(targetFeature));
    res = modPoisson(targetFeature, param, 1E-8);
end