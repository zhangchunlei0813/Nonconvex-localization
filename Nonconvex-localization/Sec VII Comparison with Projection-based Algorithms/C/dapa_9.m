clc,clear all,close all,
nf=0.0;

m=9; %anchors number
D=2; %D-dimension
L=1;

a=[-50 -59 -50 2.8 0 2.8 34 25 34;8.75 0 -8.75 2.5 0 -2.5 8.75 0 -8.75];
target=[-5;200];d=zeros(1,m);
for i=1:m
    d(i)=norm(a(:,i)-target);
end
delta=3;
kk=100000;

X_1=[];
X_2=[];
X_3=[];
X_4=[];
X_5=[];
X_6=[];
X_7=[];
X_8=[];
X_9=[];



X_difference=zeros(L,1);



for l=1:L
    tic
    r=d;
    v1=zeros(2,kk);
    w1=zeros(2,kk);
    x1=[-20 0;0 180]*ones(2,kk);
   
    v2=zeros(2,kk);
    w2=zeros(2,kk);
    x2=[-20 0;0 180]*ones(2,kk);
      
    
    v3=zeros(2,kk);
    w3=zeros(2,kk);
    x3=[-20 0;0 180]*ones(2,kk);
     
    
    v4=zeros(2,kk);
    w4=zeros(2,kk);
    x4=[-20 0;0 180]*ones(2,kk);
   
    
    v5=zeros(2,kk);
    w5=zeros(2,kk);
    x5=[-20 0;0 180]*ones(2,kk);
   
    v6=zeros(2,kk);
    w6=zeros(2,kk);
    x6=[-20 0;0 180]*ones(2,kk);
    
    v7=zeros(2,kk);
    w7=zeros(2,kk);
    x7=[-20 0;0 180]*ones(2,kk);
    
    v8=zeros(2,kk);
    w8=zeros(2,kk);
    x8=[-20 0;0 180]*ones(2,kk);
    
    v9=zeros(2,kk);
    w9=zeros(2,kk);
    x9=[-20 0;0 180]*ones(2,kk);
    
    tic
    for j=1:kk
        
        v1(:,j)=3/4*x1(:,j)+1/4*x2(:,j);
        w1(:,j)=v1(:,j)-1/(j+2)*(v1(:,j)-proj(v1(:,j),a(:,1),r(1)-delta*nf));
        x1(:,j+1)=w1(:,j)-1/(j+1)*(w1(:,j)-proj(w1(:,j),a(:,1),r(1)+delta*nf));
        %        x1(:,j+1)=v1(:,j)-1/(j+1)*(v1(:,j)-proj(v1(:,j),a(:,1),r(1+9*(l-1))+delta*nf));
        
        v2(:,j)=3/10*x2(:,j)+1/5*x5(:,j)+1/4*x1(:,j)+1/4*x3(:,j);
        w2(:,j)=v2(:,j)-1/(j+2)*(v2(:,j)-proj(v2(:,j),a(:,2),r(2)-delta*nf));
        x2(:,j+1)=w2(:,j)-1/(j+1)*(w2(:,j)-proj(w2(:,j),a(:,2),r(2)+delta*nf));
        %        x2(:,j+1)=v2(:,j)-1/(j+1)*(v2(:,j)-proj(v2(:,j),a(:,2),r(2+9*(l-1))+delta*nf));
        
        v3(:,j)=3/4*x3(:,j)+1/4*x2(:,j);
        w3(:,j)=v3(:,j)-1/(j+2)*(v3(:,j)-proj(v3(:,j),a(:,3),r(3)-delta*nf));
        x3(:,j+1)=w3(:,j)-1/(j+1)*(w3(:,j)-proj(w3(:,j),a(:,3),r(3)+delta*nf));
        %  x3(:,j+1)=v3(:,j)-1/(j+1)*(v3(:,j)-proj(v3(:,j),a(:,3),r(3+9*(l-1))+delta*nf));
        
        v4(:,j)=4/5*x4(:,j)+1/5*x5(:,j);
        w4(:,j)=v4(:,j)-1/(j+2)*(v4(:,j)-proj(v4(:,j),a(:,4),r(4)-delta*nf));
        x4(:,j+1)=w4(:,j)-1/(j+1)*(w4(:,j)-proj(w4(:,j),a(:,4),r(4)+delta*nf));
        %  x4(:,j+1)=v4(:,j)-1/(j+1)*(v4(:,j)-proj(v4(:,j),a(:,4),r(4+9*(l-1))+delta*nf));
        
        v5(:,j)=1/5*x5(:,j)+1/5*x2(:,j)+1/5*x4(:,j)+1/5*x6(:,j)+1/5*x8(:,j);
        w5(:,j)=v5(:,j)-1/(j+2)*(v5(:,j)-proj(v5(:,j),a(:,5),r(5)-delta*nf));
        x5(:,j+1)=w5(:,j)-1/(j+1)*(w5(:,j)-proj(w5(:,j),a(:,5),r(5)+delta*nf));
        % x5(:,j+1)=v5(:,j)-1/(j+1)*(v5(:,j)-proj(v5(:,j),a(:,5),r(5+9*(l-1))+delta*nf));
        
        v6(:,j)=4/5*x6(:,j)+1/5*x5(:,j);
        w6(:,j)=v6(:,j)-1/(j+2)*(v6(:,j)-proj(v6(:,j),a(:,6),r(6)-delta*nf));
        x6(:,j+1)=w6(:,j)-1/(j+1)*(w6(:,j)-proj(w6(:,j),a(:,6),r(6)+delta*nf));
        
        v7(:,j)=3/4*x7(:,j)+1/4*x8(:,j);
        w7(:,j)=v7(:,j)-1/(j+2)*(v7(:,j)-proj(v7(:,j),a(:,7),r(7)-delta*nf));
        x7(:,j+1)=w7(:,j)-1/(j+1)*(w7(:,j)-proj(w7(:,j),a(:,7),r(7)+delta*nf));
        
        v8(:,j)=3/10*x8(:,j)+1/5*x5(:,j)+1/4*x7(:,j)+1/4*x9(:,j);
        w8(:,j)=v8(:,j)-1/(j+2)*(v8(:,j)-proj(v8(:,j),a(:,8),r(8)-delta*nf));
        x8(:,j+1)=w8(:,j)-1/(j+1)*(w8(:,j)-proj(w8(:,j),a(:,8),r(8)+delta*nf));
        
        v9(:,j)=3/4*x9(:,j)+1/4*x8(:,j);
        w9(:,j)=v9(:,j)-1/(j+2)*(v9(:,j)-proj(v9(:,j),a(:,9),r(9)-delta*nf));
        x9(:,j+1)=w9(:,j)-1/(j+1)*(w9(:,j)-proj(w9(:,j),a(:,9),r(9)+delta*nf));
        
        
    end
    toc
    X_1=[X_1 x1(:,j+1)];
    X_2=[X_2 x2(:,j+1)];
    X_3=[X_3 x3(:,j+1)];
    X_4=[X_4 x4(:,j+1)];
    X_5=[X_5 x5(:,j+1)];
    X_6=[X_6 x6(:,j+1)];
    X_7=[X_7 x7(:,j+1)];
    X_8=[X_8 x8(:,j+1)];
    X_9=[X_9 x9(:,j+1)];
    
    toc
end






