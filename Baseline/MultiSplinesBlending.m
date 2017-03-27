function [offset]=MultiSplinesBlending(ftrg,trg,src,mask,posx,posy)
    src_size=size(src(:,:,1));
    S=32; % Grid Size
    new_mask=zeros(size(ftrg(:,:,1)));
    new_mask(posy:posy+src_size(1)-1, posx:posx+src_size(2)-1,:)=mask;
    grid_size=floor(size(ftrg(:,:,1))/S+2);
    %construct the Coefficents Matrix A and bx
    [A,bx]=constructMatrix(grid_size);
    %calculate outer boundary and inner boundary
    [points_outer,points_inner,g_outer,g_inner]=boundarydiffExtractor(mask,trg,src);
    points_outer(:,1)=points_outer(:,1)+posy-1;points_outer(:,2)=points_outer(:,2)+posx-1;
    points_inner(:,1)=points_inner(:,1)+posy-1;points_inner(:,2)=points_inner(:,2)+posx-1;
    %assign value to A according to the boundary points
    [A,bx]=equationConstructor(A,bx,points_outer,grid_size,g_outer,S,1);
    [A,bx]=equationConstructor(A,bx,points_inner,grid_size,g_inner,S,-1);
    %solve the equation & interpolate the membrane accoding to splines
    %control points
    for i=1:3
        x=A\bx(:,:,i);
        src_splines_points=reshape(x(1:end/2),grid_size(2),grid_size(1));
        trg_splines_points=reshape(x(end/2+1:end),grid_size(2),grid_size(1));
        src_splines_points=src_splines_points';
        trg_splines_points=trg_splines_points';
        [X,Y]=meshgrid(1:grid_size(2),1:grid_size(1));
        [Xq,Yq]=meshgrid(1:1/S:grid_size(2),1:1/S:grid_size(1));
        trg_offset=interp2(X,Y,trg_splines_points,Xq,Yq);
        src_offset=interp2(X,Y,src_splines_points,Xq,Yq);  
        offset(:,:,i)=trg_offset(1:size(ftrg,1),1:size(ftrg,2)).*~new_mask+src_offset(1:size(ftrg,1),1:size(ftrg,2)).*new_mask;    
    end
end