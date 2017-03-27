function [out] = conv2olam(a,b,mode,siz1,siz2)
%CONV2OLAM Overlap-add method of CONV2 using FFT2.
%   Y = CONV2OLAM(A,B) performs the 2-D convolution of matrices
%   A and B  using the overlap/add method and using internal parameters (FFT2
%   size and block length) which guarantee efficient execution.
%   If [ma,na] = size(A) and [mb,nb] = size(B), then size(Y) = [ma+mb-1,na+nb-1].
%
%    
%    Y = CONV2OLAM(A,B,mode,siz1,siz2) allows you to have some control over the
%    internal parameters by using a zero padding. If mode is equal to:
%      - 0 is the default value, if not specified. This value for
%          mode uses overlap/add method. Ex. conv2olam(a,b,0) is equivalent
%          to perform conv2olam(a,b).
%      - 1 the dimensions of input matrices are zero padded to their nextpow2
%          values (type 'help nextpow2' on Matlab prompt) to perform
%          FFT2 and IFFT2. No overlap is performed. Ex. conv2olam(a,b,1)
%      - 2 the dimensions of input matrices are not zero padded to perform
%          these FFT-based operations. No overlap is performed. Ex. conv2olam(a,b,2)
%      - 3 the dimensions of matrices are zero padded to fixed values
%          which must be specified to perform overlapping.
%          Ex. conv2olam(a,b,3,512,512) uses 512 X 512 matrices
%          to perform overlap/add method.
% 
%    
%    See also CONV2, FFTFILT, FFT2, IFFT2, FILTFILT.
%
%
% Please contribute if you find this software useful.
% Report bugs to luigi.rosa@tiscali.it
%
%*****************************************************************
% Luigi Rosa
% Via Centrale 27
% 67042 Civita di Bagno
% L'Aquila --- ITALY 
% email  luigi.rosa@tiscali.it
% mobile +39 340 3463208 
% http://utenti.lycos.it/matlab
%*****************************************************************
%
%

[ax,ay]=size(a);
[bx,by]=size(b);
dimx=ax+bx-1;
dimy=ay+by-1;

if (nargin<3)||(mode==0) 
    % figure out which nfftx, nffty, Lx and Ly to use
    %--------------------------------------------------------------------------    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    fftflops=zeros(13);
    fftflops(1:10,:)=[14          58         172         440        1038        2358     5264       11644       25594       55946      121644      263136     566438;...
            58         178         470        1134        2586        5738     12574       27382       59378      128274      276054      591806    1263946;...
            172         470        1170        2730        6098       13330    28858       62186      133602      286242      611498     1302394    2765458;...
            440        1134        2730        6242       13762       29794    63986      136914      292290      622658     1323346     2805490    5932322;...
            1038        2586        6098       13762       30082       64706   138210      294306      625538     1327234     2810530     5938658   12520002;...
            2358        5738       13330       29794       64706      138498   294594      624962     1323778     2799874     5911874    12458946   26203266;...
            5264       12574       28858       63986      138210      294594   624386     1320322     2788354     5881346    12386946    26044290   54659330;...
            11644       27382       62186      136914      294306      624962  1320322     2783746     5862914    12335106    25918722    54378242  113897986;...
            25594       59378      133602      292290      625538     1323778  2788354     5862914    12316674    25851906    54200834   113483266  237249538;...
            55946      128274      286242      622658     1327234     2799874  5881346    12335106    25851906    54140930   113275906   236715010  493996034];
    fftflops(11:13,:)=1.0e+009 *[0.00012164400000   0.00027605400000   0.00061149800000   0.00132334600000     0.00281053000000   0.00591187400000   0.01238694600000   0.02591872200000    0.05420083400000   0.11327590600000   0.23653990600000   0.49340621000000    1.02794445000000;...
            0.00026313600000   0.00059180600000   0.00130239400000   0.00280549000000     0.00593865800000   0.01245894600000   0.02604429000000   0.05437824200000    0.11348326600000   0.23671501000000   0.49340621000000   1.02746521800000    2.13719449800000;...
            0.00056643800000   0.00126394600000   0.00276545800000   0.00593232200000     0.01252000200000   0.02620326600000   0.05465933000000   0.11389798600000    0.23724953800000   0.49399603400000   1.02794445000000   2.13719449800000    4.43891712200000]; 
    
    
    
    
    %.....................................
    nx=1:13;
    ny=1:13;
    validsetx=find(2.^(nx-1)>bx-1);
    validsety=find(2.^(ny-1)>by-1);
    nx=nx(validsetx);
    ny=ny(validsety);
    Lx=2.^(nx-1)-bx+1;
    Ly=2.^(ny-1)-by+1;
    sizex=length(nx);
    sizey=length(ny);
    matrice=zeros(sizex,sizey);
    for ii=1:sizex
        for jj=1:sizey
            matrice(ii,jj)=ceil(dimx/Lx(ii))*ceil(dimy/Ly(jj))*fftflops(nx(ii),ny(jj));
        end
    end
    [massimo_vettore,posizione_vettore]=min(matrice);
    [massimo,posizione]=min(massimo_vettore);
    y_max=posizione;
    x_max=posizione_vettore(posizione);
    massimo;
    %.......................................
    nx2=1:13;
    ny2=1:13;
    validsetx2=find(2.^(nx2-1)>ax-1);
    validsety2=find(2.^(ny2-1)>ay-1);
    nx2=nx2(validsetx2);
    ny2=ny2(validsety2);
    Lx2=2.^(nx2-1)-ax+1;
    Ly2=2.^(ny2-1)-ay+1;
    sizex2=length(nx2);
    sizey2=length(ny2);
    matrice2=zeros(sizex2,sizey2);
    for ii=1:sizex2
        for jj=1:sizey2
            matrice2(ii,jj)=ceil(dimx/Lx2(ii))*ceil(dimy/Ly2(jj))*fftflops(nx2(ii),ny2(jj));
        end
    end
    [massimo_vettore2,posizione_vettore2]=min(matrice2);
    [massimo2,posizione2]=min(massimo_vettore2);
    y_max2=posizione2;
    x_max2=posizione_vettore2(posizione2);
    massimo2;
    %.......................................
    if massimo<massimo2
        nfftx=2^(nx(x_max)-1);
        nffty=2^(ny(y_max)-1);
        
        Lx=nfftx-bx+1;
        Ly=nffty-by+1;
        %--------------------------------------------------------------------------    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        B=fft2(b,nfftx,nffty);
        out=zeros(dimx,dimy);
        a2=a;
        a2(dimx,dimy)=0;
        
        xstart=1;
        while xstart <= dimx
            xend=min(xstart+Lx-1,dimx);
            ystart=1;
            while ystart <= dimy
                yend=min(ystart+Ly-1,dimy);
                %---------------------
                X=fft2(a2(xstart:xend,ystart:yend),nfftx,nffty);
                Y=ifft2(X.*B);
                endx=min(dimx,xstart+nfftx-1);
                endy=min(dimy,ystart+nffty-1);
                out(xstart:endx,ystart:endy)=out(xstart:endx,ystart:endy)+Y(1:(endx-xstart+1),1:(endy-ystart+1));
                ystart=ystart+Ly;
                %---------------------
            end
            xstart=xstart+Lx;
        end
        
        if ~(any(any(imag(a)))||any(any(imag(b))))
            out=real(out);
        end
        return;
    else
        nfftx=2^(nx2(x_max2)-1);
        nffty=2^(ny2(y_max2)-1);
        
        Lx2=nfftx-ax+1;
        Ly2=nffty-ay+1;
        %--------------------------------------------------------------------------    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        A=fft2(a,nfftx,nffty);
        out=zeros(dimx,dimy);
        b2=b;
        b2(dimx,dimy)=0;
        
        xstart=1;
        while xstart <= dimx
            xend=min(xstart+Lx2-1,dimx);
            ystart=1;
            while ystart <= dimy
                yend=min(ystart+Ly2-1,dimy);
                %---------------------
                X=fft2(b2(xstart:xend,ystart:yend),nfftx,nffty);
                Y=ifft2(X.*A);
                endx=min(dimx,xstart+nfftx-1);
                endy=min(dimy,ystart+nffty-1);
                out(xstart:endx,ystart:endy)=out(xstart:endx,ystart:endy)+Y(1:(endx-xstart+1),1:(endy-ystart+1));
                ystart=ystart+Ly2;
                %---------------------
            end
            xstart=xstart+Lx2;
        end
        
        if ~(any(any(imag(a)))||any(any(imag(b))))
            out=real(out);
        end
        return;
    end
else
    % mode = 1 ----> nextpow2 for both dimensions
    if mode==1
        dx=2^(nextpow2(dimx));
        dy=2^(nextpow2(dimy));
        
        out=ifft2(fft2(a,dx,dy).*fft2(b,dx,dy));
        if ~(any(any(imag(a)))||any(any(imag(b))))
            out=real(out);
        end
        out=out(1:dimx,1:dimy);
        return
    end
    
    % mode = 2 ----> the same dimensions
    if mode==2
        out=ifft2(fft2(a,dimx,dimy).*fft2(b,dimx,dimy));
        if ~(any(any(imag(a)))||any(any(imag(b))))
            out=real(out);
        end
        return
    end
    
    % mode = 3 ----> OverLap Add method with given dimensions siz1 and siz2
    if mode==3
        nfftx=siz1;
        nffty=siz2;
        
        Lx=nfftx-bx+1;
        Ly=nffty-by+1;
        
        B=fft2(b,nfftx,nffty);
        out=zeros(dimx,dimy);
        a2=a;
        a2(dimx,dimy)=0;
        
        xstart=1;
        while xstart <= dimx
            xend=min(xstart+Lx-1,dimx);
            ystart=1;
            while ystart <= dimy
                yend=min(ystart+Ly-1,dimy);
                %---------------------
                X=fft2(a2(xstart:xend,ystart:yend),nfftx,nffty);
                Y=ifft2(X.*B);
                endx=min(dimx,xstart+nfftx-1);
                endy=min(dimy,ystart+nffty-1);
                out(xstart:endx,ystart:endy)=out(xstart:endx,ystart:endy)+Y(1:(endx-xstart+1),1:(endy-ystart+1));
                ystart=ystart+Ly;
                %---------------------
            end
            xstart=xstart+Lx;
        end
        
        if ~(any(any(imag(a)))||any(any(imag(b))))
            out=real(out);
        end
        return;
        
    end
    
    
    
    
    
    
end