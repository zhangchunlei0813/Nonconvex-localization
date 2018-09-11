clc,clear all,close all,
load 'rmat_nf001.mat'
load 'X_initial_nf001.mat'
m=5;
L=1;

rho=1;
lambda=0.1;
kk2= 1000;

a1=[[5.2;1.7],[8.9;5.0],[12.1;9.0], [10.5;13.7],[ 6.3;12.2]];

target=[7.5;7.5];
d1=zeros(1,m);
for i=1:m
    d1(i)=norm(a1(:,i)-target);
end


X_difference=zeros(L,1);
N=[2 2 2 2 2];

for l=1:L
   
    a=a1;
    r=rmat(1+5*(l-1):5+5*(l-1));
    r=d1;
    for i=1:m
       x0(:,i)=X_initial(:,l);
    end
    alpha=zeros(2,m);   
    x=x0;
     tic
    for k=1:kk2
                   
        %   add alpha
        i=1;
        x(:,i)=a(:,i)+(2*r(i)+rho*norm((N(i)*x0(:,i)+x0(:,2)+x0(:,5))/2-alpha(:,i)/rho-N(i)*a(:,i)))*((N(i)*x0(:,i)+x0(:,2)+x0(:,5))/2-alpha(:,i)-N(i)*a(:,i))/((2+N(i)*rho)*norm((N(i)*x0(:,i)+x0(:,2)+x0(:,5))/2-alpha(:,i)-N(i)*a(:,i)));
        i=2;
        x(:,i)=a(:,i)+(2*r(i)+rho*norm((N(i)*x0(:,i)+x0(:,1)+x0(:,3))/2-alpha(:,i)/rho-N(i)*a(:,i)))*((N(i)*x0(:,i)+x0(:,1)+x0(:,3))/2-alpha(:,i)-N(i)*a(:,i))/((2+N(i)*rho)*norm((N(i)*x0(:,i)+x0(:,1)+x0(:,3))/2-alpha(:,i)-N(i)*a(:,i)));
        i=3;
        x(:,i)=a(:,i)+(2*r(i)+rho*norm((N(i)*x0(:,i)+x0(:,2)+x0(:,4))/2-alpha(:,i)/rho-N(i)*a(:,i)))*((N(i)*x0(:,i)+x0(:,2)+x0(:,4))/2-alpha(:,i)-N(i)*a(:,i))/((2+N(i)*rho)*norm((N(i)*x0(:,i)+x0(:,2)+x0(:,4))/2-alpha(:,i)-N(i)*a(:,i)));
        i=4;
        x(:,i)=a(:,i)+(2*r(i)+rho*norm((N(i)*x0(:,i)+x0(:,5)+x0(:,3))/2-alpha(:,i)/rho-N(i)*a(:,i)))*((N(i)*x0(:,i)+x0(:,5)+x0(:,3))/2-alpha(:,i)-N(i)*a(:,i))/((2+N(i)*rho)*norm((N(i)*x0(:,i)+x0(:,5)+x0(:,3))/2-alpha(:,i)-N(i)*a(:,i)));
        i=5;
        x(:,i)=a(:,i)+(2*r(i)+rho*norm((N(i)*x0(:,i)+x0(:,1)+x0(:,4))/2-alpha(:,i)/rho-N(i)*a(:,i)))*((N(i)*x0(:,i)+x0(:,1)+x0(:,4))/2-alpha(:,i)-N(i)*a(:,i))/((2+N(i)*rho)*norm((N(i)*x0(:,i)+x0(:,1)+x0(:,4))/2-alpha(:,i)-N(i)*a(:,i)));
        
       
        x0=x;
        
        alpha(:,1)=alpha(:,1)+rho*(N(1)*x(:,1)-x(:,2)-x(:,5))/2;
        alpha(:,2)=alpha(:,2)+rho*(N(2)*x(:,2)-x(:,1)-x(:,3))/2;
        alpha(:,3)=alpha(:,3)+rho*(N(3)*x(:,3)-x(:,2)-x(:,4))/2;
        alpha(:,4)=alpha(:,4)+rho*(N(4)*x(:,4)-x(:,5)-x(:,3))/2;
        alpha(:,5)=alpha(:,5)+rho*(N(5)*x(:,5)-x(:,4)-x(:,1))/2;
 
    end
    
    toc
end
% INC=sqrt(sum(X_difference)/L);
% ERR1=sqrt(X_sum1/L);
% ERR2=sqrt(X_sum2/L);
% ERR3=sqrt(X_sum3/L);
% ERR4=sqrt(X_sum4/L);
% ERR5=sqrt(X_sum5/L);
% 
% XX1=zeros(1,199);
% XX2=zeros(1,199);
% XX3=zeros(1,199);
% XX4=zeros(1,199);
% XX5=zeros(1,199);
% R1=zeros(1,199);
% R2=zeros(1,199);
% R3=zeros(1,199);
% R4=zeros(1,199);
% R5=zeros(1,199);
% XX_average=zeros(1,199);
% for t=1:199
%     XX1(:,t)=norm(X1(:,t+1)-X1(:,t));
%     R1(:,t)=norm(X1(:,t)-a(:,1));
% end
% for t=1:199
%     XX2(:,t)=norm(X2(:,t+1)-X2(:,t));
%     R2(:,t)=norm(X2(:,t)-a(:,2));
% end
% for t=1:199
%     XX3(:,t)=norm(X3(:,t+1)-X3(:,t));
%     R3(:,t)=norm(X3(:,t)-a(:,3));
% end
% for t=1:199
%     XX4(:,t)=norm(X4(:,t+1)-X4(:,t));
%     R4(:,t)=norm(X4(:,t)-a(:,4));
% end
% for t=1:199
%     XX5(:,t)=norm(X5(:,t+1)-X5(:,t));
%     R5(:,t)=norm(X5(:,t)-a(:,5));
% end

% figure
% plot(target(1,1),target(2,1),'k+','linewidth',2,'markersize',10)
% hold on
% for l=1:L
%     plot(X_1(1,l), X_1(2,l), 'r*','linewidth',1,'markersize',5)
%     hold on
% end
% for l=1:L
%     plot(X_2(1,l), X_2(2,l), 'bo','linewidth',1,'markersize',5)
%     hold on
% end
% for l=1:L
%     plot(X_3(1,l), X_3(2,l), 'gs','linewidth',1,'markersize',5)
%     hold on
% end
% for l=1:L
%     plot(X_4(1,l), X_4(2,l), 'yd','linewidth',1,'markersize',5)
%     hold on
% end
% for l=1:L
%     plot(X_5(1,l), X_5(2,l), 'cv','linewidth',1,'markersize',5)
%     hold on
% end
%
% axis([0 15 0 15])
% xlabel('X')
% ylabel('Y')
% grid on

