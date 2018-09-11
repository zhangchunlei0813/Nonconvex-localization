clc,clear all,close all
m=5;
L=1000; % number of trials
kk=1000; % iteration 

% collect the solution obtained at each trial
X_1=[];
X_2=[];
X_3=[];
X_4=[];
X_5=[];

N=[2 3 2 3 2];
RHO=zeros(5*L,kk+1);

for l=1:L
    rho=1*ones(m,kk+1); % initialize rho
    a=[[5.2;1.7],[8.9;5.0],[12.1;9.0], [10.5;13.7],[ 6.3;12.2]]; 
    % position of sensors
    for i=1:m
        t= 2*pi*rand(1,1);
        b=15*rand(1,1);
        a(:,i)=b*[sin(t);cos(t)];
    end
    %true target position
    target=[7.5;7.5];
    d=zeros(1,m);
    for i=1:m
        d(i)=norm(a(:,i)-target);
    end
    r=d;
    tic
    x0=zeros(2,m);
    theta=pi/2;
    for i=1:m
        x0(:,i)=[a(1,i)+r(i)*cos(theta) a(2,i)+r(i)*sin(theta)];
    end
    lambda=zeros(2,m);
    
    % record the iteration of each sensor's estimate value 
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
        
        % update lambda
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
        
        
        % update x
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
        
        % update rho
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
    
    
   
    X_1=[X_1 x(:,1)];    
    X_2=[X_2 x(:,2)];
    X_3=[X_3 x(:,3)];
    X_4=[X_4 x(:,4)];
    X_5=[X_5 x(:,5)];
    toc
    
    % record the iteration of rho
    RHO(5*(l-1)+1:5*l,:)=rho;
end






