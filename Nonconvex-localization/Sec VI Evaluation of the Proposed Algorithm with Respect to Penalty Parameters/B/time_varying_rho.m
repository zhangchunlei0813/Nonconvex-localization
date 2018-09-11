clc,clear all,close all
load 'rmat_nf001.mat'
load 'X_initial_nf001.mat'
m=5;
L=1000;
X_sum1=0;
X_sum2=0;
X_sum3=0;
X_sum4=0;
X_sum5=0;
kk=1000;

X_difference=zeros(L,1);
N=[2 3 2 3 2];
X_1=[];
X_2=[];
X_3=[];
X_4=[];
X_5=[];

RHO=zeros(5*L,kk+1);

for l=1:L
    rho=1*ones(m,kk+1);
    a=[[5.2;1.7],[8.9;5.0],[12.1;9.0], [10.5;13.7],[ 6.3;12.2]];
    target=[7.5;7.5];
    r=rmat(1+5*(l-1):5+5*(l-1));
    tic
    for i=1:m
        x0(:,i)=X_initial(:,l);
    end

    lambda=zeros(2,m);
    
    X1=[];
    X2=[];
    X3=[];
    X4=[];
    X5=[];
    
    x=x0;
    
    for k=1:kk
        X1=[X1 x0(:,1)];
        X2=[X2 x0(:,2)];
        X3=[X3 x0(:,3)];
        X4=[X4 x0(:,4)];
        X5=[X5 x0(:,5)];
        if k==1
            rho0=rho(1:m,k);
        else
            rho0=rho(1:m,k-1);
        end
        
        %update lambda
        i=1;
        lambda(:,1)=lambda(:,1)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5));
        
        i=2;
        lambda(:,2)=lambda(:,2)+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3))+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        
        i=3;
        lambda(:,3)=lambda(:,3)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        
        i=4;
        lambda(:,4)=lambda(:,4)+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3));
        
        i=5;
        lambda(:,5)=lambda(:,5)+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4))+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1));%         end
        
        
        %update x
        i=1;
        mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))-N(i)*rho(i,k)*a(:,i)-lambda(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);

        i=2;
        mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1))+1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-lambda(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        
        i=3;
        mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-lambda(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        
        i=4;
        mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3))+1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))-N(i)*rho(i,k)*a(:,i)-lambda(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        
        i=5;
        mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-lambda(:,i);
        x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        
        %update rho
        vv=10^-12;
        epsilon=0.99;
        i=1;
        if norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5)),inf)<=vv
            rho(1,k+1)=max([rho(1,k),rho(2,k),rho(5,k)]);
        else
            if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5)),inf)
                rho(1,k+1)=max([rho(1,k),rho(2,k),rho(5,k)]);
            else
                rho(1,k+1)=1.0001*max([rho(1,k),rho(2,k),rho(5,k)]);
            end
        end
        i=2;
        if norm(rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=vv;
            rho(2,k+1)=max([rho(1,k),rho(2,k),rho(3,k),rho(4,k)]);
        else
            if norm(rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1))+rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3))+rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
                rho(2,k+1)=max([rho(1,k),rho(2,k),rho(3,k),rho(4,k)]);
                
            else
                rho(2,k+1)=1.0001*max([rho(1,k),rho(2,k),rho(3,k),rho(4,k)]);
            end
        end
        i=3;
        if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=vv
            rho(3,k+1)=max([rho(2,k),rho(3,k),rho(4,k)]);
        else
            if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
                rho(3,k+1)=max([rho(2,k),rho(3,k),rho(4,k)]);
            else
                rho(3,k+1)=1.0001*max([rho(2,k),rho(3,k),rho(4,k)]);
            end
        end
        
        i=4;
        if norm(rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3)),inf)<=vv
            rho(4,k+1)=max([rho(2,k),rho(3,k),rho(4,k),rho(5,k)]);
        else
            if norm(rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3)),inf)<=epsilon*norm(rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5))+rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3)),inf)
                rho(4,k+1)=max([rho(2,k),rho(3,k),rho(4,k),rho(5,k)]);
            else
                rho(4,k+1)=1.0001*max([rho(2,k),rho(3,k),rho(4,k),rho(5,k)]);
            end
        end
        
        i=5;
        if norm(rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4))+rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1)),inf)<=vv
            rho(5,k+1)=max([rho(1,k),rho(4,k),rho(5,k)]);
        else
            if norm(rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4))+rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1)),inf)<=epsilon*norm(rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4))+rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1)),inf)
                rho(5,k+1)=max([rho(1,k),rho(4,k),rho(5,k)]);
                
            else
                rho(5,k+1)=1.0001*max([rho(1,k),rho(4,k),rho(5,k)]);
                
            end
        end
        x0=x;
        
    end
    
    
    X_sum1=X_sum1+(norm(x(:,1)-target))^2;
    X_1=[X_1 x(:,1)];
    X_sum2=X_sum2+(norm(x(:,2)-target))^2;
    X_2=[X_2 x(:,2)];
    X_sum3=X_sum3+(norm(x(:,3)-target))^2;
    X_3=[X_3 x(:,3)];
    X_sum4=X_sum4+(norm(x(:,4)-target))^2;
    X_4=[X_4 x(:,4)];
    X_sum5=X_sum5+(norm(x(:,5)-target))^2;
    X_5=[X_5 x(:,5)];
    
    for i=1:m-1
        for j=i+1:m
            X_difference(l)=X_difference(l)+(norm(x(:,i)-x(:,j)))^2;
        end
    end
    toc
    RHO(5*(l-1)+1:5*l,:)=rho;
end
INC=sqrt(sum(X_difference)/L);
ERR1=sqrt(X_sum1/L);
ERR2=sqrt(X_sum2/L);
ERR3=sqrt(X_sum3/L);
ERR4=sqrt(X_sum4/L);
ERR5=sqrt(X_sum5/L);

figure
plot(target(1,1),target(2,1),'k+','linewidth',2,'markersize',10)
hold on
for l=1:L
    plot(X_1(1,l), X_1(2,l), 'r*','linewidth',1,'markersize',5)
    hold on
end
for l=1:L
    plot(X_2(1,l), X_2(2,l), 'bo','linewidth',1,'markersize',5)
    hold on
end
for l=1:L
    plot(X_3(1,l), X_3(2,l), 'gs','linewidth',1,'markersize',5)
    hold on
end
for l=1:L
    plot(X_4(1,l), X_4(2,l), 'yd','linewidth',1,'markersize',5)
    hold on
end
for l=1:L
    plot(X_5(1,l), X_5(2,l), 'cv','linewidth',1,'markersize',5)
    hold on
end

axis([0 15 0 15])
xlabel('X')
ylabel('Y')
grid on

