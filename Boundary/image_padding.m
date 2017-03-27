
function [padding_src,padding_trg,padding_mask_src,padding_mask_trg]=image_padding(src,trg,mask)
    padding_size=power(2,ceil(log2(size(src(:,:,1)))));
    padding_src=zeros(padding_size(1),padding_size(2),3);
    padding_trg=zeros(padding_size(1),padding_size(2),3);
    padding_src(1:size(src(:,:,1),1),1:size(src(:,:,1),2),:)=src;
    padding_trg(1:size(src(:,:,1),1),1:size(src(:,:,1),2),:)=trg;
    mask_src = repmat(mask,[1,1,3]);
    mask_trg = 1-mask_src;
    padding_mask_src=zeros(padding_size(1),padding_size(2),3);
    padding_mask_trg=ones(padding_size(1),padding_size(2),3);
    padding_mask_src(1:size(src(:,:,1),1),1:size(src(:,:,1),2),:)=mask_src;
    padding_mask_trg(1:size(src(:,:,1),1),1:size(src(:,:,1),2),:)=mask_trg;    
end