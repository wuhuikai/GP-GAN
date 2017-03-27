function [naive_img]=NaiveBlending(trg,src,mask)
    naive_img=trg;
    naive_img(repmat(mask,[1,1,3]))=src(repmat(mask,[1,1,3]));
end