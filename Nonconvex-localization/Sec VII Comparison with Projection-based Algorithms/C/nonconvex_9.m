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
m=9;
X_sum1=0;
X_sum2=0;
X_sum3=0;
X_sum4=0;
X_sum5=0;
X_sum6=0;
X_sum7=0;
X_sum8=0;
X_sum9=0;

a=[[-58.93;0],[-50;8.75],[-50;-8.75],...
    [0;0],[2.8;2.5],[2.8;-2.5],...
    [25;0],[33.93;8.75],[33.93;-8.75]];

lambda=0.1;
nf=0.0;
kk2= 500;
L=10000;
X_difference=zeros(L,1);
epsilon=0.99;

target=[-5;200];
rtrue=zeros(m,1);
N=[3 1 1 4 1 1 3 1 1];
X_1=[];
X_2=[];
X_3=[];
X_4=[];
X_5=[];
X_6=[];
X_7=[];
X_8=[];
X_9=[];
HX=[];
x0=zeros(2,m);
rho=ones(m,kk2+1);
for i=1:m
    rtrue(i)=norm(a(:,i)-target);
end
for l=1:L
    tic

    x0(:,1)=X_initial1(:,l);
    x0(:,2)=X_initial2(:,l);
    x0(:,3)=X_initial3(:,l);
    x0(:,4)=X_initial4(:,l);
    x0(:,5)=X_initial5(:,l);
    x0(:,6)=X_initial6(:,l);
    x0(:,7)=X_initial7(:,l);
    x0(:,8)=X_initial8(:,l);
    x0(:,9)=X_initial9(:,l);
    
    alpha=zeros(2,m);

    r=rtrue+nf.*randn(m,1);
    x=x0;
    
    for k=1:kk2
        if k==1
            rho0=rho(1:m,k);
        else
            rho0=rho(1:m,k-1);
        end
        
        i=1;
        mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        i=2;
        mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        i=3;
        mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        i=4;
        mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1))+1/(rho0(i)+rho0(6))*(rho0(i)*x(:,i)+rho0(6)*x(:,6))+1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5))+1/(rho0(i)+rho0(7))*(rho0(i)*x(:,i)+rho0(7)*x(:,7)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        i=5;
        mu=rho(i,k)*(1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
         i=6;
        mu=rho(i,k)*(1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
         i=7;
        mu=rho(i,k)*(1/(rho0(i)+rho0(8))*(rho0(i)*x(:,i)+rho0(8)*x(:,8))+1/(rho0(i)+rho0(9))*(rho0(i)*x(:,i)+rho0(9)*x(:,9))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
         i=8;
        mu=rho(i,k)*(1/(rho0(i)+rho0(7))*(rho0(i)*x(:,i)+rho0(7)*x(:,7)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
         i=9;
        mu=rho(i,k)*(1/(rho0(i)+rho0(7))*(rho0(i)*x(:,i)+rho0(7)*x(:,7)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);

        i=1;
        alpha(:,1)=alpha(:,1)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3))+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        i=2;
        alpha(:,2)=alpha(:,2)+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1));
        i=3;
        alpha(:,3)=alpha(:,3)+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1));
        i=4;
        alpha(:,4)=alpha(:,4)+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(i,k)*rho(6,k)/(rho(i,k)+rho(6,k))*(x(:,i)-x(:,6))+rho(i,k)*rho(7,k)/(rho(i,k)+rho(7,k))*(x(:,i)-x(:,7))+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1));
        i=5;
        alpha(:,5)=alpha(:,5)+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        i=6;
        alpha(:,6)=alpha(:,6)+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        i=7;
        alpha(:,7)=alpha(:,7)+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4))+rho(i,k)*rho(8,k)/(rho(i,k)+rho(8,k))*(x(:,i)-x(:,8))+rho(i,k)*rho(9,k)/(rho(i,k)+rho(9,k))*(x(:,i)-x(:,9));
        i=8;
        alpha(:,8)=alpha(:,8)+rho(i,k)*rho(7,k)/(rho(i,k)+rho(7,k))*(x(:,i)-x(:,7));
        i=9;
        alpha(:,9)=alpha(:,9)+rho(i,k)*rho(7,k)/(rho(i,k)+rho(7,k))*(x(:,i)-x(:,7));
        
        
        i=1;
        if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3))+rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
            rho(1,k+1)=max([rho(1,k),rho(2,k),rho(3,k),rho(4,k)]);
        else
            rho(1,k+1)=1.0001*max([rho(1,k),rho(2,k),rho(2,k),rho(4,k)]);
        end
        i=2;
        if norm(rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1)),inf)<=epsilon*norm(rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1)),inf)
            rho(2,k+1)=max([rho(1,k),rho(2,k)]);
        else
            rho(2,k+1)=1.0001*max([rho(1,k),rho(2,k)]);
        end
        i=3;
        if norm(rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1)),inf)<=epsilon*norm(rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1)),inf)
            rho(3,k+1)=max([rho(1,k),rho(3,k)]);
        else
            rho(3,k+1)=1.0001*max([rho(1,k),rho(3,k)]);
        end
        i=4;
        if norm(rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(6,k)/(rho(i,k)+rho(6,k))*(x(:,i)-x(:,6))+rho(7,k)/(rho(i,k)+rho(7,k))*(x(:,i)-x(:,7))+rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1)),inf)<=epsilon*norm(rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5))+rho0(6)/(rho0(i)+rho0(6))*(x0(:,i)-x0(:,6))+rho0(7)/(rho0(i)+rho0(7))*(x0(:,i)-x0(:,7))+rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1)),inf)
            rho(4,k+1)=max([rho(1,k),rho(5,k),rho(6,k),rho(7,k),rho(4,k)]);
        else
            rho(4,k+1)=1.0001*max([rho(1,k),rho(5,k),rho(6,k),rho(7,k),rho(4,k)]);
        end
        i=5;
        if norm(rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
            rho(5,k+1)=max([rho(4,k),rho(5,k)]);
        else
            rho(5,k+1)=1.0001*max([rho(4,k),rho(5,k)]);
        end
        i=6;
        if norm(rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
            rho(6,k+1)=max([rho(4,k),rho(6,k)]);
        else
            rho(6,k+1)=1.0001*max([rho(4,k),rho(6,k)]);
        end
        i=7;
        if norm(rho(8,k)/(rho(i,k)+rho(8,k))*(x(:,i)-x(:,8))+rho(9,k)/(rho(i,k)+rho(9,k))*(x(:,i)-x(:,9))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(8)/(rho0(i)+rho0(8))*(x0(:,i)-x0(:,8))+rho0(9)/(rho0(i)+rho0(9))*(x0(:,i)-x0(:,9))+rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
            rho(7,k+1)=max([rho(7,k),rho(8,k),rho(9,k),rho(4,k)]);
        else
            rho(7,k+1)=1.0001*max([rho(7,k),rho(8,k),rho(9,k),rho(4,k)]);
        end
        i=8;
        if norm(rho(7,k)/(rho(i,k)+rho(7,k))*(x(:,i)-x(:,7)),inf)<=epsilon*norm(rho0(7)/(rho0(i)+rho0(7))*(x0(:,i)-x0(:,7)),inf)
            rho(8,k+1)=max([rho(7,k),rho(8,k)]);
        else
            rho(8,k+1)=1.0001*max([rho(7,k),rho(8,k)]);
        end
        i=9;
        if norm(rho(7,k)/(rho(i,k)+rho(7,k))*(x(:,i)-x(:,7)),inf)<=epsilon*norm(rho0(7)/(rho0(i)+rho0(7))*(x0(:,i)-x0(:,7)),inf)
            rho(9,k+1)=max([rho(7,k),rho(9,k)]);
        else
            rho(9,k+1)=1.0001*max([rho(7,k),rho(9,k)]);
        end
        
        x0=x;
        
    end
    
    

    X_1=[X_1 x(:,1)];
    X_2=[X_2 x(:,2)];
    X_3=[X_3 x(:,3)];
    X_4=[X_4 x(:,4)];
    X_5=[X_5 x(:,5)];
    X_6=[X_6 x(:,6)];
    X_7=[X_7 x(:,7)];
    X_8=[X_8 x(:,8)];
    X_9=[X_9 x(:,9)];
    
    toc
end

AA= length(find(X_1(2,:)<0));
