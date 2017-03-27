function [A,bx] = equationConstructor(A,bx,points,trg_grid,g,S,flag)
    km=floor(points/S);
    target_coefficient=[km(:,1)*trg_grid(2)+km(:,2)+1,km(:,1)*trg_grid(2)+km(:,2)+1+1, ...
        (km(:,1)+1)*trg_grid(2)+km(:,2)+1,(km(:,1)+1)*trg_grid(2)+km(:,2)+1+1];
    bsplines=[bspline_basis((points(:,1)-km(:,1)*S)/S).*bspline_basis((points(:,2)-km(:,2)*S)/S), ...
        bspline_basis((points(:,1)-km(:,1)*S)/S).*bspline_basis((points(:,2)-(km(:,2)+1)*S)/S), ...
        bspline_basis((points(:,1)-(km(:,1)+1)*S)/S).*bspline_basis((points(:,2)-(km(:,2))*S)/S), ...
        bspline_basis((points(:,1)-(km(:,1)+1)*S)/S).*bspline_basis((points(:,2)-(km(:,2)+1)*S)/S)
        ]; 
    source_coefficient=trg_grid(1)*trg_grid(2)+[km(:,1)*trg_grid(2)+km(:,2)+1,km(:,1)*trg_grid(2)+km(:,2)+1+1, ...
        (km(:,1)+1)*trg_grid(2)+km(:,2)+1,(km(:,1)+1)*trg_grid(2)+km(:,2)+1+1];
    idx=zeros(0,2);
    adt=zeros(0,1);
    for i=1:4
        x=repmat(source_coefficient(:,i),4,1);
        y=reshape(source_coefficient,[],1);
        z=reshape(target_coefficient,[],1);
        idx=[idx;[x,y];[x,z]];
        x_=repmat(bsplines(:,i),4,1);
        y_=-1*reshape(bsplines,[],1);
        z_=-y_;
        adt=[adt;x_.*y_;x_.*z_];
        for j=1:size(source_coefficient(:,i),1)
            for channel=1:3
                bx(source_coefficient(j,i),1,channel)=bx(source_coefficient(j,i),1,channel)+flag*bsplines(j,i)*g(j,1,channel); 
            end
        end
    end
      for i=1:4
          x=repmat(target_coefficient(:,i),4,1);
          y=reshape(source_coefficient,[],1);
          z=reshape(target_coefficient,[],1);
          idx=[idx;[x,y];[x,z]];
         for j=1:size(source_coefficient(:,i),1)
             for channel=1:3
                bx(target_coefficient(j,i),1,channel)=bx(target_coefficient(j,i),1,channel)-flag*bsplines(j,i)*g(j,1,channel); 
             end
        end
          %bx(target_coefficient(:,i))=bx(target_coefficient(:,i))-flag*bsplines(:,i).*g;
      end
    adt=[adt;-adt];
    idx=sub2ind(size(A),idx(:,1),idx(:,2));
    for i=1:size(adt,1)
        A(idx(i))=A(idx(i))+adt(i);
    end
end

