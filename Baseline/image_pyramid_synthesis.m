function [pyr]=image_pyramid_synthesis(pyr_src,pyr_trg,mask_src,mask_trg,level)
    for p = 1:level
        mask_src = imresize(mask_src,[size(pyr_src{p},1),size(pyr_src{p},2)]);
        mask_trg = imresize(mask_trg,[size(pyr_src{p},1),size(pyr_src{p},2)]);
        pyr_synthesis{p} = pyr_src{p}.*mask_src + pyr_trg{p}.*mask_trg;
    end
    for p=length(pyr_synthesis)-1:-1:1
        pyr_synthesis{p}=pyr_synthesis{p}+impyramid_(pyr_synthesis{p+1},'expand');
    end
    pyr=pyr_synthesis{1};
end