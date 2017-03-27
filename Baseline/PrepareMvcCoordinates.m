function [distance,vertex,p,t,angle_left,angle_right]=PrepareMvcCoordinates(chi,mask)
        %get counter-clockwise vertexs 
        [vertex_x,vertex_y]=find(chi);
        vertex=[vertex_x,vertex_y];
        [~,vertex_order] = sort([angle(complex(vertex_x-mean(vertex_x),vertex_y-mean(vertex_y)))]);
        vertex=vertex(vertex_order,:);
        %[coordinates_x,coordinates_y]=find(mask);
        %coordinates=[coordinates_x,coordinates_y];
        idx=cell2mat(bwboundaries(mask));        
        [p,t]=mesh2d((idx));
       % boundarydiff_vector=boundarydiff(sub2ind(size(error),vertex(:,1),vertex(:,2)));     
        distance=pdist2(p,vertex);
        
        vertex_horizon=reshape(vertex',1,[]);
        angle_vector(:,:,1)=vertex_horizon;
        angle_vector(:,:,2)=circshift(vertex_horizon,[0,2])';
        angle_vector(:,:,3)=circshift(vertex_horizon,[0,-2])';

        vector_p=repmat(p,[1,size(vertex,1)]);
        vector_boundary=repmat(angle_vector,[size(p,1),1,1]);
        vector_middle=vector_boundary(:,:,1)-vector_p;
        vector_left=vector_boundary(:,:,2)-vector_p;
        vector_right=vector_boundary(:,:,3)-vector_p;
        
        angle_left=vector_middle.*vector_left;
        angle_right=vector_middle.*vector_right;

        norm_angle_left=vector_left.^2;
        norm_angle_right=vector_right.^2;
        norm_angle_middle=vector_middle.^2;
        
        angle_left=angle_left(:,1:2:end)+angle_left(:,2:2:end);
        angle_right=angle_right(:,1:2:end)+angle_right(:,2:2:end);
        norm_angle_left=norm_angle_left(:,1:2:end)+norm_angle_left(:,2:2:end);
        norm_angle_right=norm_angle_right(:,1:2:end)+norm_angle_right(:,2:2:end);
        norm_angle_middle=norm_angle_middle(:,1:2:end)+norm_angle_middle(:,2:2:end);

        angle_left=real(acos(angle_left./((norm_angle_left.*norm_angle_middle).^(1/2))));
        angle_right=real(acos(angle_right./((norm_angle_right.*norm_angle_middle).^(1/2))));
        
end