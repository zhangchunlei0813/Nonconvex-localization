clc,clear all,close all,
load 'rmat_nf001.mat'
load 'X_initial_nf001.mat'
nf=0.01;

m=5; %anchors number
D=2; %D-dimension
L=1000;

a=[5.2 8.9 12.1 10.5 6.3;1.7 5.0 9.0 13.7 12.2];
target=[7.5;7.5];d=zeros(1,m);
for i=1:m
    d(i)=norm(a(:,i)-target);
end
delta=3;
kk=1000;


tic
for l=1:L
    r=rmat(1+5*(l-1):5+5*(l-1));
    v1=zeros(2,kk);
    w1=zeros(2,kk);
    x1=[-10 0;0 180]*ones(2,kk);
    x1(:,1)=X_initial(:,l);
    
    v2=zeros(2,kk);
    w2=zeros(2,kk);
    x2=[-10 0;0 180]*ones(2,kk);
    x2(:,1)=X_initial(:,l);
    
    
    v3=zeros(2,kk);
    w3=zeros(2,kk);
    x3=[-10 0;0 180]*ones(2,kk);
    x3(:,1)=X_initial(:,l);
    
    
    v4=zeros(2,kk);
    w4=zeros(2,kk);
    x4=[-10 0;0 180]*ones(2,kk);
    x4(:,1)=X_initial(:,l);
    
    
    v5=zeros(2,kk);
    w5=zeros(2,kk);
    x5=[-10 0;0 180]*ones(2,kk);
    x5(:,1)=X_initial(:,l);
    
    tic
    for j=1:kk 
        v1(:,j)=1/3*x1(:,j)+1/3*x2(:,j)+1/3*x5(:,j);
        w1(:,j)=v1(:,j)-1/(j+2)*(v1(:,j)-proj(v1(:,j),a(:,1),r(1)-delta*nf));
        x1(:,j+1)=w1(:,j)-1/(j+1)*(w1(:,j)-proj(w1(:,j),a(:,1),r(1)+delta*nf));
        %        x1(:,j+1)=v1(:,j)-1/(j+1)*(v1(:,j)-proj(v1(:,j),a(:,1),r(1+9*(l-1))+delta*nf));
        
        v2(:,j)=1/3*x2(:,j)+1/3*x1(:,j)+1/3*x3(:,j);
        w2(:,j)=v2(:,j)-1/(j+2)*(v2(:,j)-proj(v2(:,j),a(:,2),r(2)-delta*nf));
        x2(:,j+1)=w2(:,j)-1/(j+1)*(w2(:,j)-proj(w2(:,j),a(:,2),r(2)+delta*nf));
        %        x2(:,j+1)=v2(:,j)-1/(j+1)*(v2(:,j)-proj(v2(:,j),a(:,2),r(2+9*(l-1))+delta*nf));
        
        v3(:,j)=1/3*x3(:,j)+1/3*x2(:,j)+1/3*x4(:,j);
        w3(:,j)=v3(:,j)-1/(j+2)*(v3(:,j)-proj(v3(:,j),a(:,3),r(3)-delta*nf));
        x3(:,j+1)=w3(:,j)-1/(j+1)*(w3(:,j)-proj(w3(:,j),a(:,3),r(3)+delta*nf));
        %  x3(:,j+1)=v3(:,j)-1/(j+1)*(v3(:,j)-proj(v3(:,j),a(:,3),r(3+9*(l-1))+delta*nf));
        
        v4(:,j)=1/3*x4(:,j)+1/3*x5(:,j)+1/3*x3(:,j);
        w4(:,j)=v4(:,j)-1/(j+2)*(v4(:,j)-proj(v4(:,j),a(:,4),r(4)-delta*nf));
        x4(:,j+1)=w4(:,j)-1/(j+1)*(w4(:,j)-proj(w4(:,j),a(:,4),r(4)+delta*nf));
        %  x4(:,j+1)=v4(:,j)-1/(j+1)*(v4(:,j)-proj(v4(:,j),a(:,4),r(4+9*(l-1))+delta*nf));
        
        v5(:,j)=1/3*x5(:,j)+1/3*x4(:,j)+1/3*x1(:,j);
        w5(:,j)=v5(:,j)-1/(j+2)*(v5(:,j)-proj(v5(:,j),a(:,5),r(5)-delta*nf));
        x5(:,j+1)=w5(:,j)-1/(j+1)*(w5(:,j)-proj(w5(:,j),a(:,5),r(5)+delta*nf));
        % x5(:,j+1)=v5(:,j)-1/(j+1)*(v5(:,j)-proj(v5(:,j),a(:,5),r(5+9*(l-1))+delta*nf));
        
    end
    toc
end
 toc
