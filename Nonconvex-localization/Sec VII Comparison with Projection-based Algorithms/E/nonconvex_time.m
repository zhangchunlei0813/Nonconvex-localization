clc,clear all,close all
load 'rmat_nf001.mat'
load 'X_initial_nf001.mat'
m=5;
L=1000;

kk2= 1000;

a=[[5.2;1.7],[8.9;5.0],[12.1;9.0], [10.5;13.7],[ 6.3;12.2]];

target=[7.5;7.5];
d=zeros(1,m);
for i=1:m
    d(i)=norm(a(:,i)-target);
end
X_difference=zeros(L,1);
N=[2 2 2 2 2];


RHO=zeros(5*L,kk2+1);
tic
for l=1:L
    rho=1*ones(m,kk2+1);
    r=rmat(1+5*(l-1):5+5*(l-1));
    tic
    for i=1:m
        x0(:,i)=X_initial(:,l);
    end
    
    alpha=zeros(2,m);
    
    X1=[];
    X2=[];
    X3=[];
    X4=[];
    X5=[];
    
    x=x0;
    tic
    for k=1:kk2
        r=rmat(1+5*(l-1):5+5*(l-1));
        if k==1
            rho0=rho(1:m,k);
        else
            rho0=rho(1:m,k-1);
        end
        
        
        %   add alpha
        i=1;
        mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        %       hx(:,i)=2*(x(:,i)-a(:,i))-2*r(i)*(x(:,i)-a(:,i))/norm(x(:,i)-a(:,i))+ rho(i,k)*N(i)*x(:,i)-rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))+alpha(:,i);
        %       HX=[HX hx(:,1)];
        i=2;
        mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1))+1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        i=3;
        mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        i=4;
        mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(3)*x(:,3))+1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        i=5;
        mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        
        i=1;
        alpha(:,1)=alpha(:,1)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5));
        
        
        i=2;
        alpha(:,2)=alpha(:,2)+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3));
        
        i=3;
        alpha(:,3)=alpha(:,3)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        
        i=4;
        alpha(:,4)=alpha(:,4)+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3));
        
        i=5;
        alpha(:,5)=alpha(:,5)+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4))+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1));
        
        epsilon=0.99;
        i=1;
        
        if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5)),inf)
            rho(1,k+1)=max([rho(1,k),rho(2,k),rho(5,k)]);
            
        else
            rho(1,k+1)=1.0001*max([rho(1,k),rho(2,k),rho(5,k)]);
            
        end
        i=2;
        
        if norm(rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3)),inf)<=epsilon*norm(rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1))+rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3)),inf)
            rho(2,k+1)=max([rho(1,k),rho(2,k),rho(3,k)]);
            
        else
            rho(2,k+1)=1.0001*max([rho(1,k),rho(2,k),rho(3,k)]);
            
        end
        i=3;
        
        if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
            rho(3,k+1)=max([rho(2,k),rho(3,k),rho(4,k)]);
            
        else
            rho(3,k+1)=1.0001*max([rho(2,k),rho(3,k),rho(4,k)]);
            
        end
        
        i=4;
        
        if norm(rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3)),inf)<=epsilon*norm(rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5))+rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3)),inf)
            rho(4,k+1)=max([rho(3,k),rho(4,k),rho(5,k)]);
            
        else
            rho(4,k+1)=1.0001*max([rho(3,k),rho(4,k),rho(5,k)]);
            
        end
        
        i=5;
        if norm(rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4))+rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1)),inf)<=epsilon*norm(rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4))+rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1)),inf)
            rho(5,k+1)=max([rho(1,k),rho(4,k),rho(5,k)]);
            
        else
            rho(5,k+1)=1.0001*max([rho(1,k),rho(4,k),rho(5,k)]);
            
        end
        
        x0=x;
        RHO(5*(l-1)+1:5*l,:)=rho;
    end
    toc
end
toc
%INC=sqrt(sum(X_difference)/L);
% ERR1=sqrt(X_sum1/L);
% ERR2=sqrt(X_sum2/L);
% ERR3=sqrt(X_sum3/L);
% ERR4=sqrt(X_sum4/L);
% ERR5=sqrt(X_sum5/L);

