clc,clear all,close all,
load 'rmat_nf001.mat'
load 'X_initial_nf001.mat'
m=5; %anchors number
D=2; %D-dimension
L=1000;

a=[5.2 8.9 12.1 10.5 6.3;1.7 5.0 9.0 13.7 12.2];
target=[7.5;7.5];
d1=zeros(1,m);
for i=1:m
    d1(i)=norm(a(:,i)-target);
end

kk=1000;

BETA=[];
for j=1:L
    tic
     r=rmat(1+5*(j-1):5+5*(j-1));
    % r(5)=1.7119;
    alpha(:,1)=X_initial(:,j);
    for k=1:kk
        beta=alpha(:,k);
        for i=1:m
          alpha(:,k+1)=a(:,i)+r(i)*(beta-a(:,i))/norm(beta-a(:,i));
            beta=alpha(:,k+1);
        end
    end
  toc
end

