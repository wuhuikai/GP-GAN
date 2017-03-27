function [res]=MVCBlending(trg,src,mask)
    %extract boundary & calculate boundary difference
    h = fspecial('laplacian', 0);
    chi = imfilter(double(mask),h);
    chi(chi<0) = 0;
    chi(chi>0) = 1;  
    erf=trg - src;
    res_copy=trg.*~repmat(mask,[1,1,3])+src.*repmat(mask,[1,1,3]);
    boundarydiff=erf.*repmat(chi,[1,1,3]);
    %pre-calculation for mvc 
    [distance,vertex,p,t,angle_left,angle_right]=PrepareMvcCoordinates(chi,mask);
    W=(tan(angle_left/2)+tan(angle_right/2))./distance;
    Lamda=W./repmat(sum(W,2),[1,size(W,2)]);
    vector_bdiff=zeros(size(W,1),size(W,2),3);
    for i=1:3
        bdiff=boundarydiff(:,:,i);
        bdiff=bdiff(sub2ind(size(erf),vertex(:,1),vertex(:,2)));
        bdiff=repmat(bdiff',[size(Lamda,1),1]);
        vector_bdiff(:,:,i)=bdiff;
    end
    vector_offset=sum(vector_bdiff.*repmat(Lamda,[1,1,3]),2);
    [pi(:,1),pi(:,2)]=find(ones(size(mask,1),size(mask,2)));
    i=1:size(mask,1)*size(mask,2);i=i';
    DT=delaunayTriangulation(p);
    for i=1:3
        fi=tri_interp(p,t,vector_offset(:,:,i),pi,i,DT);
        offset_mvc(:,:,i)=reshape(fi,size(mask)); 
    end
    res=res_copy+offset_mvc.*repmat(mask,[1,1,3]);
end