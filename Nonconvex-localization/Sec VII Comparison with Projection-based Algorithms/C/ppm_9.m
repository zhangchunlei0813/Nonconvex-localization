clc,clear all,close all,
load 'X_initial1_9_nf0_500_10000.mat'
load 'X_initial2_9_nf0_500_10000.mat'
load 'X_initial3_9_nf0_500_10000.mat'
load 'X_initial4_9_nf0_500_10000.mat'
load 'X_initial5_9_nf0_500_10000.mat'
load 'X_initial6_9_nf0_500_10000.mat'
load 'X_initial7_9_nf0_500_10000.mat'
load 'X_initial8_9_nf0_500_10000.mat'
load 'X_initial9_9_nf0_500_10000.mat'
m=9; %anchors number
D=2; %D-dimension
L=10000;

a=[-58.93 -50 -50 0 2.8 2.8 25 33.93 33.93; 0 8.75 -8.75 0 2.5 -2.5 0 8.75 -8.75];
target=[-5;200];
r=zeros(1,m);

for i=1:m
    r(i)=norm(a(:,i)-target);
end

kk=500;
Xmat=[];


for l=1:L
    tic
    x0=[-50 0;100 0]*ones(2,kk);

    x0(:,1)=X_initial1(:,l);
    x0(:,2)=X_initial2(:,l);
    x0(:,3)=X_initial3(:,l);
    x0(:,4)=X_initial4(:,l);
    x0(:,5)=X_initial5(:,l);
    x0(:,6)=X_initial6(:,l);
    x0(:,7)=X_initial7(:,l);
    x0(:,8)=X_initial8(:,l);
    x0(:,9)=X_initial9(:,l); 
    

%     
    for j=1:kk
        x=0;
        for i=1:m
            x=x+a(:,i)+r(i)*(x0(:,j)-a(:,i))/norm(x0(:,j)-a(:,i));
        end
        x=x/m;
        x0(:,j+1)=x;
        
    end
    toc
    Xmat=[Xmat x];
end

Xmat_sum=0;
for nn=1:L
    Xmat_sum=Xmat_sum+(norm(Xmat(:,nn)-target))^2;
end
Xmat_average=sqrt(Xmat_sum/L);
AA= length(find(Xmat(2,:)<0));




