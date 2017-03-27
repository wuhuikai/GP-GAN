function [pyr_src,pyr_trg]=image_pyramid_analysis(padding_src,padding_trg,level)
    pyr_src{1}=padding_src;
    for p=2:level
        pyr_src{p}=impyramid(pyr_src{p-1},'reduce');
    end
    for p = 1:level-1
        pyr_src{p} = pyr_src{p}-impyramid_(pyr_src{p+1},'expand');
    end
    pyr_trg{1}=padding_trg;
    for p=2:level
        pyr_trg{p}=impyramid(pyr_trg{p-1},'reduce');
    end
    for p = 1:level-1
        pyr_trg{p} = pyr_trg{p}-impyramid_(pyr_trg{p+1},'expand');
    end
end