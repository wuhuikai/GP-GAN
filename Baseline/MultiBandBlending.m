function [res]=MultiBandBlending(trg,src,mask,level)
    [padding_src,padding_trg,padding_mask_src,padding_mask_trg]=image_padding(src,trg,mask);
    [pyr_src,pyr_trg]=image_pyramid_analysis(padding_src,padding_trg,level);
    h = fspecial('gauss',5,3); 
    mask_src = imfilter(padding_mask_src,h,'replicate');
    mask_trg = imfilter(padding_mask_trg,h,'replicate');
    res=image_pyramid_synthesis(pyr_src,pyr_trg,mask_src,mask_trg,level);
    res=res(1:size(src(:,:,1),1),1:size(src(:,:,1),2),:);
end