
function [points_return_o,points_return_i,go,gi]=boundarydiffExtractor(mask,trg,src)
     h=fspecial('laplacian',0);
     outer_chi=imfilter(double(mask),h);
     outer_chi(outer_chi<0)=0;outer_chi(outer_chi>0)=1;
     inner_chi=zeros(size(src(:,:,1)));
     idx=cell2mat(bwboundaries(mask));
     inner_chi(sub2ind(size(src),idx(:,1),idx(:,2)))=1;
     
     error=trg-src;
     
     left_x=outer_chi.*circshift(inner_chi,[0,-1]);left_right_x=circshift(left_x,[0,1]);
     right_x=inner_chi.*circshift(outer_chi,[0,-1]);right_left_x=circshift(right_x,[0,1]);
     up_y=outer_chi.*circshift(inner_chi,[-1,0]);up_down_y=circshift(up_y,[1,0]);
     down_y=inner_chi.*circshift(outer_chi,[-1,0]);down_up_y=circshift(down_y,[1,0]);
     

     points_o=[find(left_x);find(up_y);];
     points_i=[find(right_x);find(down_y)];
     
     for i=1:3
        error_channel=error(:,:,i);
        go(:,:,i)=[error_channel(find(left_x))+error_channel(find(left_right_x)); ...
         error_channel(find(up_y))+error_channel(find(up_down_y))]; ...
        gi(:,:,i)= [-error_channel(find(right_x))-error_channel(find(right_left_x)); ...
         -error_channel(find(down_y))-error_channel(find(down_up_y))];    
     end
     go=go/2;gi=gi/2;
     [points_return_o(:,1),points_return_o(:,2)]=ind2sub(size(src(:,:,1)),points_o);
     [points_return_i(:,1),points_return_i(:,2)]=ind2sub(size(src(:,:,1)),points_i);
end