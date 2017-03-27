function [A,bx]=constructMatrix(trg_grid)
    %construct the source image splines part
    t=trg_grid(1)*trg_grid(2);
    A_trg=eye(t);
    A_trg=-4*A_trg+triu(circshift(A_trg,[0,trg_grid(2)]),0)+tril(circshift(A_trg,[0,-trg_grid(2)]),0) ...
        +triu(circshift(A_trg,[0,1]),0)+tril(circshift(A_trg,[0,-1]),0);
    zero_term=trg_grid(2):trg_grid(2):trg_grid(1)*trg_grid(2)-1;
    zero_term=[[zero_term'+1;zero_term'],[zero_term';zero_term'+1]];
    A_trg(sub2ind(size(A_trg),zero_term(:,1),zero_term(:,2)))=0;
    %composite the source and target  to A & construct bx
    A=zeros(2*t);
    A(1:t,1:t)=A_trg;
    A(t+1:2*t,t+1:2*t)=A_trg;  
    bx=zeros(2*t,1,3);
end