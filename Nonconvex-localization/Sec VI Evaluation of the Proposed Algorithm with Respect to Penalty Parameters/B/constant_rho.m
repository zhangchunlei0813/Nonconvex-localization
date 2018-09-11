clc,clear all,close all,
load 'rmat_nf2.mat'
load 'X_initial_nf2.mat'
m=5;
L=1000;
X_sum1=0;
X_sum2=0;
X_sum3=0;
X_sum4=0;
X_sum5=0;

rho=1;
kk= 1000;

a=[[5.2;1.7],[8.9;5.0],[12.1;9.0], [10.5;13.7],[ 6.3;12.2]];
target=[7.5;7.5];


X_difference=zeros(L,1);
N=[2 3 2 3 2];
X_1=[];
X_2=[];
X_3=[];
X_4=[];
X_5=[];

for l=1:L
    tic
    r=rmat(1+5*(l-1):5+5*(l-1));
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
        
        %update x
        i=1;
        x(:,i)=a(:,i)+(2*r(i)+rho*norm((N(i)*x0(:,i)+x0(:,2)+x0(:,5))/2-lambda(:,i)/rho-N(i)*a(:,i)))*((N(i)*x0(:,i)+x0(:,2)+x0(:,5))/2-lambda(:,i)-N(i)*a(:,i))/((2+N(i)*rho)*norm((N(i)*x0(:,i)+x0(:,2)+x0(:,5))/2-lambda(:,i)-N(i)*a(:,i)));
        i=2;
        x(:,i)=a(:,i)+(2*r(i)+rho*norm((N(i)*x0(:,i)+x0(:,1)+x0(:,3)+x0(:,4))/2-lambda(:,i)/rho-N(i)*a(:,i)))*((N(i)*x0(:,i)+x0(:,1)+x0(:,3)+x0(:,4))/2-lambda(:,i)-N(i)*a(:,i))/((2+N(i)*rho)*norm((N(i)*x0(:,i)+x0(:,1)+x0(:,3)+x0(:,4))/2-lambda(:,i)-N(i)*a(:,i)));
        i=3;
        x(:,i)=a(:,i)+(2*r(i)+rho*norm((N(i)*x0(:,i)+x0(:,2)+x0(:,4))/2-lambda(:,i)/rho-N(i)*a(:,i)))*((N(i)*x0(:,i)+x0(:,2)+x0(:,4))/2-lambda(:,i)-N(i)*a(:,i))/((2+N(i)*rho)*norm((N(i)*x0(:,i)+x0(:,2)+x0(:,4))/2-lambda(:,i)-N(i)*a(:,i)));
        i=4;
        x(:,i)=a(:,i)+(2*r(i)+rho*norm((N(i)*x0(:,i)+x0(:,2)+x0(:,5)+x0(:,3))/2-lambda(:,i)/rho-N(i)*a(:,i)))*((N(i)*x0(:,i)+x0(:,2)+x0(:,5)+x0(:,3))/2-lambda(:,i)-N(i)*a(:,i))/((2+N(i)*rho)*norm((N(i)*x0(:,i)+x0(:,2)+x0(:,5)+x0(:,3))/2-lambda(:,i)-N(i)*a(:,i)));
        i=5;
        x(:,i)=a(:,i)+(2*r(i)+rho*norm((N(i)*x0(:,i)+x0(:,1)+x0(:,4))/2-lambda(:,i)/rho-N(i)*a(:,i)))*((N(i)*x0(:,i)+x0(:,1)+x0(:,4))/2-lambda(:,i)-N(i)*a(:,i))/((2+N(i)*rho)*norm((N(i)*x0(:,i)+x0(:,1)+x0(:,4))/2-lambda(:,i)-N(i)*a(:,i)));
        
        x0=x;
        
        lambda(:,1)=lambda(:,1)+rho*(N(1)*x(:,1)-x(:,2)-x(:,5))/2;
        lambda(:,2)=lambda(:,2)+rho*(N(2)*x(:,2)-x(:,1)-x(:,3)-x(:,4))/2;
        lambda(:,3)=lambda(:,3)+rho*(N(3)*x(:,3)-x(:,2)-x(:,4))/2;
        lambda(:,4)=lambda(:,4)+rho*(N(4)*x(:,4)-x(:,2)-x(:,5)-x(:,3))/2;
        lambda(:,5)=lambda(:,5)+rho*(N(5)*x(:,5)-x(:,4)-x(:,1))/2;
 
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

