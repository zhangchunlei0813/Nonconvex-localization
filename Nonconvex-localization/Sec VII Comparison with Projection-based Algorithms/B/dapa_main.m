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



X_difference=zeros(L,1);
X_sum1=0;
X_sum2=0;
X_sum3=0;
X_sum4=0;
X_sum5=0;


X_1=[];
X_2=[];
X_3=[];
X_4=[];
X_5=[];

Xmat1_1=[];
Xmat1_10=[];
Xmat1_20=[];
Xmat1_30=[];
Xmat1_40=[];
Xmat1_50=[];
Xmat1_100=[];
Xmat1_500=[];


Xmat2_1=[];
Xmat2_10=[];
Xmat2_20=[];
Xmat2_30=[];
Xmat2_40=[];
Xmat2_50=[];
Xmat2_100=[];
Xmat2_500=[];


Xmat3_1=[];
Xmat3_10=[];
Xmat3_20=[];
Xmat3_30=[];
Xmat3_40=[];
Xmat3_50=[];
Xmat3_100=[];
Xmat3_500=[];


Xmat4_1=[];
Xmat4_10=[];
Xmat4_20=[];
Xmat4_30=[];
Xmat4_40=[];
Xmat4_50=[];
Xmat4_100=[];
Xmat4_500=[];


Xmat5_1=[];
Xmat5_10=[];
Xmat5_20=[];
Xmat5_30=[];
Xmat5_40=[];
Xmat5_50=[];
Xmat5_100=[];
Xmat5_500=[];



for l=1:L
    tic
    r=rmat(1+5*(l-1):5+5*(l-1));
    v1=zeros(2,kk);
    w1=zeros(2,kk);
    x1(:,1)=X_initial(:,l);
    
    v2=zeros(2,kk);
    w2=zeros(2,kk);
    x2(:,1)=X_initial(:,l);
        
    v3=zeros(2,kk);
    w3=zeros(2,kk);
    x3(:,1)=X_initial(:,l);
        
    v4=zeros(2,kk);
    w4=zeros(2,kk);
    x4(:,1)=X_initial(:,l);
        
    v5=zeros(2,kk);
    w5=zeros(2,kk);
    x5(:,1)=X_initial(:,l);
    
    tic
    for j=1:kk
        
        v1(:,j)=5/12*x1(:,j)+1/4*x2(:,j)+1/3*x5(:,j);
        w1(:,j)=v1(:,j)-1/(j+2)*(v1(:,j)-proj(v1(:,j),a(:,1),r(1)-delta*nf));
        x1(:,j+1)=w1(:,j)-1/(j+1)*(w1(:,j)-proj(w1(:,j),a(:,1),r(1)+delta*nf));
        %        x1(:,j+1)=v1(:,j)-1/(j+1)*(v1(:,j)-proj(v1(:,j),a(:,1),r(1+9*(l-1))+delta*nf));
        
        v2(:,j)=1/4*x2(:,j)+1/4*x1(:,j)+1/4*x3(:,j)+1/4*x4(:,j);
        w2(:,j)=v2(:,j)-1/(j+2)*(v2(:,j)-proj(v2(:,j),a(:,2),r(2)-delta*nf));
        x2(:,j+1)=w2(:,j)-1/(j+1)*(w2(:,j)-proj(w2(:,j),a(:,2),r(2)+delta*nf));
        %        x2(:,j+1)=v2(:,j)-1/(j+1)*(v2(:,j)-proj(v2(:,j),a(:,2),r(2+9*(l-1))+delta*nf));
        
        v3(:,j)=1/2*x3(:,j)+1/4*x2(:,j)+1/4*x4(:,j);
        w3(:,j)=v3(:,j)-1/(j+2)*(v3(:,j)-proj(v3(:,j),a(:,3),r(3)-delta*nf));
        x3(:,j+1)=w3(:,j)-1/(j+1)*(w3(:,j)-proj(w3(:,j),a(:,3),r(3)+delta*nf));
        %  x3(:,j+1)=v3(:,j)-1/(j+1)*(v3(:,j)-proj(v3(:,j),a(:,3),r(3+9*(l-1))+delta*nf));
        
        v4(:,j)=1/4*x4(:,j)+1/4*x2(:,j)+1/4*x5(:,j)+1/4*x3(:,j);
        w4(:,j)=v4(:,j)-1/(j+2)*(v4(:,j)-proj(v4(:,j),a(:,4),r(4)-delta*nf));
        x4(:,j+1)=w4(:,j)-1/(j+1)*(w4(:,j)-proj(w4(:,j),a(:,4),r(4)+delta*nf));
        %  x4(:,j+1)=v4(:,j)-1/(j+1)*(v4(:,j)-proj(v4(:,j),a(:,4),r(4+9*(l-1))+delta*nf));
        
        v5(:,j)=5/12*x5(:,j)+1/4*x4(:,j)+1/3*x1(:,j);
        w5(:,j)=v5(:,j)-1/(j+2)*(v5(:,j)-proj(v5(:,j),a(:,5),r(5)-delta*nf));
        x5(:,j+1)=w5(:,j)-1/(j+1)*(w5(:,j)-proj(w5(:,j),a(:,5),r(5)+delta*nf));
        % x5(:,j+1)=v5(:,j)-1/(j+1)*(v5(:,j)-proj(v5(:,j),a(:,5),r(5+9*(l-1))+delta*nf));
        
        if j==100
            Xmat1_100=[Xmat1_100 x1(:,j+1)];
            Xmat2_100=[Xmat2_100 x2(:,j+1)];
            Xmat3_100=[Xmat3_100 x3(:,j+1)];
            Xmat4_100=[Xmat4_100 x4(:,j+1)];
            Xmat5_100=[Xmat5_100 x5(:,j+1)];
        end
        if j==10
            Xmat1_10=[Xmat1_10 x1(:,j+1)];
            Xmat2_10=[Xmat2_10 x2(:,j+1)];
            Xmat3_10=[Xmat3_10 x3(:,j+1)];
            Xmat4_10=[Xmat4_10 x4(:,j+1)];
            Xmat5_10=[Xmat5_10 x5(:,j+1)];
        end
        if j==500
            Xmat1_500=[Xmat1_500 x1(:,j+1)];
            Xmat2_500=[Xmat2_500 x2(:,j+1)];
            Xmat3_500=[Xmat3_500 x3(:,j+1)];
            Xmat4_500=[Xmat4_500 x4(:,j+1)];
            Xmat5_500=[Xmat5_500 x5(:,j+1)];
        end
        if j==50
            Xmat1_50=[Xmat1_50 x1(:,j+1)];
            Xmat2_50=[Xmat2_50 x2(:,j+1)];
            Xmat3_50=[Xmat3_50 x3(:,j+1)];
            Xmat4_50=[Xmat4_50 x4(:,j+1)];
            Xmat5_50=[Xmat5_50 x5(:,j+1)];
        end
        if j==1
            Xmat1_1=[Xmat1_1 x1(:,j+1)];
            Xmat2_1=[Xmat2_1 x2(:,j+1)];
            Xmat3_1=[Xmat3_1 x3(:,j+1)];
            Xmat4_1=[Xmat4_1 x4(:,j+1)];
            Xmat5_1=[Xmat5_1 x5(:,j+1)];
        end
        if j==20
            Xmat1_20=[Xmat1_20 x1(:,j+1)];
            Xmat2_20=[Xmat2_20 x2(:,j+1)];
            Xmat3_20=[Xmat3_20 x3(:,j+1)];
            Xmat4_20=[Xmat4_20 x4(:,j+1)];
            Xmat5_20=[Xmat5_20 x5(:,j+1)];
        end
        if j==30
            Xmat1_30=[Xmat1_30 x1(:,j+1)];
            Xmat2_30=[Xmat2_30 x2(:,j+1)];
            Xmat3_30=[Xmat3_30 x3(:,j+1)];
            Xmat4_30=[Xmat4_30 x4(:,j+1)];
            Xmat5_30=[Xmat5_30 x5(:,j+1)];
        end
        if j==40
            Xmat1_40=[Xmat1_40 x1(:,j+1)];
            Xmat2_40=[Xmat2_40 x2(:,j+1)];
            Xmat3_40=[Xmat3_40 x3(:,j+1)];
            Xmat4_40=[Xmat4_40 x4(:,j+1)];
            Xmat5_40=[Xmat5_40 x5(:,j+1)];
        end
        
    end
    toc
    X_1=[X_1 x1(:,j+1)];
    X_2=[X_2 x2(:,j+1)];
    X_3=[X_3 x3(:,j+1)];
    X_4=[X_4 x4(:,j+1)];
    X_5=[X_5 x5(:,j+1)];
    
    
    X_sum1=X_sum1+(norm(x1(:,j+1)-target))^2;
    X_sum2=X_sum2+(norm(x2(:,j+1)-target))^2;
    X_sum3=X_sum3+(norm(x3(:,j+1)-target))^2;
    X_sum4=X_sum4+(norm(x4(:,j+1)-target))^2;
    X_sum5=X_sum5+(norm(x5(:,j+1)-target))^2;
    
    toc
end

ERR1=sqrt(X_sum1/L);
ERR2=sqrt(X_sum2/L);
ERR3=sqrt(X_sum3/L);
ERR4=sqrt(X_sum4/L);
ERR5=sqrt(X_sum5/L);



Xmat_sum1_10=0;
for nn=1:L
    Xmat_sum1_10=Xmat_sum1_10+(norm(Xmat1_10(:,nn)-target))^2;
end
Xmat_average1_10=sqrt(Xmat_sum1_10/L);

Xmat_sum1_1=0;
for nn=1:L
    Xmat_sum1_1=Xmat_sum1_1+(norm(Xmat1_1(:,nn)-target))^2;
end
Xmat_average1_1=sqrt(Xmat_sum1_1/L);

Xmat_sum1_100=0;
for nn=1:L
    Xmat_sum1_100=Xmat_sum1_100+(norm(Xmat1_100(:,nn)-target))^2;
end
Xmat_average1_100=sqrt(Xmat_sum1_100/L);

Xmat_sum1_500=0;
for nn=1:L
    Xmat_sum1_500=Xmat_sum1_500+(norm(Xmat1_500(:,nn)-target))^2;
end
Xmat_average1_500=sqrt(Xmat_sum1_500/L);

Xmat_sum1_50=0;
for nn=1:L
    Xmat_sum1_50=Xmat_sum1_50+(norm(Xmat1_50(:,nn)-target))^2;
end
Xmat_average1_50=sqrt(Xmat_sum1_50/L);

Xmat_sum1_20=0;
for nn=1:L
    Xmat_sum1_20=Xmat_sum1_20+(norm(Xmat1_20(:,nn)-target))^2;
end
Xmat_average1_20=sqrt(Xmat_sum1_20/L);

Xmat_sum1_30=0;
for nn=1:L
    Xmat_sum1_30=Xmat_sum1_30+(norm(Xmat1_30(:,nn)-target))^2;
end
Xmat_average1_30=sqrt(Xmat_sum1_30/L);

Xmat_sum1_40=0;
for nn=1:L
    Xmat_sum1_40=Xmat_sum1_40+(norm(Xmat1_40(:,nn)-target))^2;
end
Xmat_average1_40=sqrt(Xmat_sum1_40/L);

%.....................................................................
Xmat_sum2_10=0;
for nn=1:L
    Xmat_sum2_10=Xmat_sum2_10+(norm(Xmat2_10(:,nn)-target))^2;
end
Xmat_average2_10=sqrt(Xmat_sum2_10/L);

Xmat_sum2_1=0;
for nn=1:L
    Xmat_sum2_1=Xmat_sum2_1+(norm(Xmat2_1(:,nn)-target))^2;
end
Xmat_average2_1=sqrt(Xmat_sum2_1/L);

Xmat_sum2_100=0;
for nn=1:L
    Xmat_sum2_100=Xmat_sum2_100+(norm(Xmat2_100(:,nn)-target))^2;
end
Xmat_average2_100=sqrt(Xmat_sum2_100/L);

Xmat_sum2_500=0;
for nn=1:L
    Xmat_sum2_500=Xmat_sum2_500+(norm(Xmat2_500(:,nn)-target))^2;
end
Xmat_average2_500=sqrt(Xmat_sum2_500/L);

Xmat_sum2_50=0;
for nn=1:L
    Xmat_sum2_50=Xmat_sum2_50+(norm(Xmat2_50(:,nn)-target))^2;
end
Xmat_average2_50=sqrt(Xmat_sum2_50/L);

Xmat_sum2_20=0;
for nn=1:L
    Xmat_sum2_20=Xmat_sum2_20+(norm(Xmat2_20(:,nn)-target))^2;
end
Xmat_average2_20=sqrt(Xmat_sum2_20/L);

Xmat_sum2_30=0;
for nn=1:L
    Xmat_sum2_30=Xmat_sum2_30+(norm(Xmat2_30(:,nn)-target))^2;
end
Xmat_average2_30=sqrt(Xmat_sum2_30/L);

Xmat_sum2_40=0;
for nn=1:L
    Xmat_sum2_40=Xmat_sum2_40+(norm(Xmat2_40(:,nn)-target))^2;
end
Xmat_average2_40=sqrt(Xmat_sum2_40/L);

%.....................................................................
Xmat_sum3_10=0;
for nn=1:L
    Xmat_sum3_10=Xmat_sum3_10+(norm(Xmat3_10(:,nn)-target))^2;
end
Xmat_average3_10=sqrt(Xmat_sum3_10/L);

Xmat_sum3_1=0;
for nn=1:L
    Xmat_sum3_1=Xmat_sum3_1+(norm(Xmat3_1(:,nn)-target))^2;
end
Xmat_average3_1=sqrt(Xmat_sum3_1/L);

Xmat_sum3_100=0;
for nn=1:L
    Xmat_sum3_100=Xmat_sum3_100+(norm(Xmat3_100(:,nn)-target))^2;
end
Xmat_average3_100=sqrt(Xmat_sum3_100/L);

Xmat_sum3_500=0;
for nn=1:L
    Xmat_sum3_500=Xmat_sum3_500+(norm(Xmat3_500(:,nn)-target))^2;
end
Xmat_average3_500=sqrt(Xmat_sum3_500/L);

Xmat_sum3_50=0;
for nn=1:L
    Xmat_sum3_50=Xmat_sum3_50+(norm(Xmat3_50(:,nn)-target))^2;
end
Xmat_average3_50=sqrt(Xmat_sum3_50/L);

Xmat_sum3_20=0;
for nn=1:L
    Xmat_sum3_20=Xmat_sum3_20+(norm(Xmat3_20(:,nn)-target))^2;
end
Xmat_average3_20=sqrt(Xmat_sum3_20/L);

Xmat_sum3_30=0;
for nn=1:L
    Xmat_sum3_30=Xmat_sum3_30+(norm(Xmat3_30(:,nn)-target))^2;
end
Xmat_average3_30=sqrt(Xmat_sum3_30/L);

Xmat_sum3_40=0;
for nn=1:L
    Xmat_sum3_40=Xmat_sum3_40+(norm(Xmat3_40(:,nn)-target))^2;
end
Xmat_average3_40=sqrt(Xmat_sum3_40/L);

%.....................................................................
Xmat_sum4_10=0;
for nn=1:L
    Xmat_sum4_10=Xmat_sum4_10+(norm(Xmat4_10(:,nn)-target))^2;
end
Xmat_average4_10=sqrt(Xmat_sum4_10/L);

Xmat_sum4_1=0;
for nn=1:L
    Xmat_sum4_1=Xmat_sum4_1+(norm(Xmat4_1(:,nn)-target))^2;
end
Xmat_average4_1=sqrt(Xmat_sum4_1/L);

Xmat_sum4_100=0;
for nn=1:L
    Xmat_sum4_100=Xmat_sum4_100+(norm(Xmat4_100(:,nn)-target))^2;
end
Xmat_average4_100=sqrt(Xmat_sum4_100/L);

Xmat_sum4_500=0;
for nn=1:L
    Xmat_sum4_500=Xmat_sum4_500+(norm(Xmat4_500(:,nn)-target))^2;
end
Xmat_average4_500=sqrt(Xmat_sum4_500/L);

Xmat_sum4_50=0;
for nn=1:L
    Xmat_sum4_50=Xmat_sum4_50+(norm(Xmat4_50(:,nn)-target))^2;
end
Xmat_average4_50=sqrt(Xmat_sum4_50/L);

Xmat_sum4_20=0;
for nn=1:L
    Xmat_sum4_20=Xmat_sum4_20+(norm(Xmat4_20(:,nn)-target))^2;
end
Xmat_average4_20=sqrt(Xmat_sum4_20/L);

Xmat_sum4_30=0;
for nn=1:L
    Xmat_sum4_30=Xmat_sum4_30+(norm(Xmat4_30(:,nn)-target))^2;
end
Xmat_average4_30=sqrt(Xmat_sum4_30/L);

Xmat_sum4_40=0;
for nn=1:L
    Xmat_sum4_40=Xmat_sum4_40+(norm(Xmat4_40(:,nn)-target))^2;
end
Xmat_average4_40=sqrt(Xmat_sum4_40/L);

%.....................................................................
Xmat_sum5_10=0;
for nn=1:L
    Xmat_sum5_10=Xmat_sum5_10+(norm(Xmat5_10(:,nn)-target))^2;
end
Xmat_average5_10=sqrt(Xmat_sum5_10/L);

Xmat_sum5_1=0;
for nn=1:L
    Xmat_sum5_1=Xmat_sum5_1+(norm(Xmat5_1(:,nn)-target))^2;
end
Xmat_average5_1=sqrt(Xmat_sum5_1/L);

Xmat_sum5_100=0;
for nn=1:L
    Xmat_sum5_100=Xmat_sum5_100+(norm(Xmat5_100(:,nn)-target))^2;
end
Xmat_average5_100=sqrt(Xmat_sum5_100/L);

Xmat_sum5_500=0;
for nn=1:L
    Xmat_sum5_500=Xmat_sum5_500+(norm(Xmat5_500(:,nn)-target))^2;
end
Xmat_average5_500=sqrt(Xmat_sum5_500/L);

Xmat_sum5_50=0;
for nn=1:L
    Xmat_sum5_50=Xmat_sum5_50+(norm(Xmat5_50(:,nn)-target))^2;
end
Xmat_average5_50=sqrt(Xmat_sum5_50/L);

Xmat_sum5_20=0;
for nn=1:L
    Xmat_sum5_20=Xmat_sum5_20+(norm(Xmat5_20(:,nn)-target))^2;
end
Xmat_average5_20=sqrt(Xmat_sum5_20/L);

Xmat_sum5_30=0;
for nn=1:L
    Xmat_sum5_30=Xmat_sum5_30+(norm(Xmat5_30(:,nn)-target))^2;
end
Xmat_average5_30=sqrt(Xmat_sum5_30/L);

Xmat_sum5_40=0;
for nn=1:L
    Xmat_sum5_40=Xmat_sum5_40+(norm(Xmat5_40(:,nn)-target))^2;
end
Xmat_average5_40=sqrt(Xmat_sum5_40/L);





