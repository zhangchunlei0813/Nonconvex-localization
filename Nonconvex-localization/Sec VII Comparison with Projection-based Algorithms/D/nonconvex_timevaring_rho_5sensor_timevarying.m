clc,clear all,close all,
load 'rmat_nf001.mat'
load 'X_initial_nf001.mat'
load 'NN_05.mat'
load 'TT_05.mat'
m=5;
L=1000;
X_sum1=0;
X_sum2=0;
X_sum3=0;
X_sum4=0;
X_sum5=0;

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

lambda=0.1;
nf=0;
kk2= 1000;
kk1=300;
a=[[5.2;1.7],[8.9;5.0],[12.1;9.0], [10.5;13.7],[ 6.3;12.2]];
target=[7.5;7.5];
d=zeros(1,m);
for i=1:m
    d(i)=norm(a(:,i)-target);
end

N=[0 0 0 0 0];

X_difference=zeros(L,1);
X_1=[];
X_2=[];
X_3=[];
X_4=[];
X_5=[];

HX=[];


for l=1:L
    rho=ones(m,kk2+1);
    r=rmat(1+5*(l-1):5+5*(l-1));
    tic
    x0=[-50 0;100 0]*ones(2,m);
    for i=1:m
        x0(:,i)=X_initial(:,l);
    end
    alpha=zeros(2,m);
    X1=[];
    X2=[];
    X3=[];
    X4=[];
    X5=[];
    
    x1=[];
    x=x0;
    
    for k=1:kk2
        %          t=rand(1,1);
        %         if t<0.5
        %             t1=0;
        %         else
        %             t1=1;
        %         end
        %
        %         t=rand(1,1);
        %         if t<0.5
        %             t2=0;
        %         else
        %             t2=1;
        %         end
        %
        %         t=rand(1,1);
        %         if t<0.5
        %             t3=0;
        %         else
        %             t3=1;
        %         end
        %
        %         t=rand(1,1);
        %         if t<0.5
        %             t4=0;
        %         else
        %             t4=1;
        %         end
        %
        %         t=rand(1,1);
        %         if t<0.5
        %             t5=0;
        %         else
        %             t5=1;
        %         end
        %
        %
        %         t=rand(1,1);
        %         if t<0.5
        %             t6=0;
        %         else
        %             t6=1;
        %         end
        %
        %
        %         N(1)=length(find([t1 t5]>0.5));
        %         N(2)=length(find([t1 t2 t6]>0.5));
        %         N(3)=length(find([t2 t3]>0.5));
        %         N(4)=length(find([t3 t4 t6]>0.5));
        %         N(5)=length(find([t4 t5]>0.5));
        N=NN(k,1+m*(l-1):m*l);
        T=TT(k,1+6*(l-1):6*l);
        t1=T(1);
        t2=T(2);
        t3=T(3);
        t4=T(4);
        t5=T(5);
        t6=T(6);
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
        
        
        
        i=1;
        if t1==1&&t5==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t1==1&&t5==0
            mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if  t1==0&&t5==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if  t1==0&&t5==0
            x(:,i)=x(:,i);
        end
        
        i=2;
        if t1==1&&t2==1&&t6==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1))+1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t1==1&&t2==1&&t6==0
            mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1))+1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t1==1&&t2==0&&t6==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t1==0&&t2==1&&t6==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t1==1&&t2==0&&t6==0
            mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t1==0&&t2==1&&t6==0
            mu=rho(i,k)*(1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t1==0&&t2==0&&t6==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t1==0&&t2==0&&t6==0
            x(:,i)=x(:,i);
        end
        
        
        i=3;
        if t2==1&&t3==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t2==1&&t3==0
            mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t2==0&&t3==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t2==0&&t3==0
            x(:,i)=x(:,i);
        end
        
        i=4;
        if t3==1&&t4==1&&t6==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3))+1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t3==1&&t4==1&&t6==0
            mu=rho(i,k)*(1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3))+1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t3==1&&t4==0&&t6==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t3==0&&t4==1&&t6==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2))+1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t3==1&&t4==0&&t6==0
            mu=rho(i,k)*(1/(rho0(i)+rho0(3))*(rho0(i)*x(:,i)+rho0(3)*x(:,3)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t3==0&&t4==1&&t6==0
            mu=rho(i,k)*(1/(rho0(i)+rho0(5))*(rho0(i)*x(:,i)+rho0(5)*x(:,5)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t3==0&&t4==0&&t6==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(2))*(rho0(i)*x(:,i)+rho0(2)*x(:,2)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t3==0&&t4==0&&t6==0
            x(:,i)=x(:,i);
        end
        
        i=5;
        if t4==1&&t5==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1))+1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t4==1&&t5==0
            mu=rho(i,k)*(1/(rho0(i)+rho0(4))*(rho0(i)*x(:,i)+rho0(4)*x(:,4)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t4==0&&t5==1
            mu=rho(i,k)*(1/(rho0(i)+rho0(1))*(rho0(i)*x(:,i)+rho0(1)*x(:,1)))-N(i)*rho(i,k)*a(:,i)-alpha(:,i);
            x(:,i)=a(:,i)+(2*r(i)+norm(mu))/(2+N(i)*rho(i,k))*mu/norm(mu);
        end
        if t4==0&&t5==0
            x(:,i)=x(:,i);
        end
        
        
        i=1;
        if t1==1&&t5==1
            alpha(:,1)=alpha(:,1)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5));
        end
        if t1==1&&t5==0
            alpha(:,1)=alpha(:,1)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2));
        end
        if t1==0&&t5==1
            alpha(:,1)=alpha(:,1)+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5));
        end
        if t1==0&&t5==0
            alpha(:,1)=alpha(:,1);
        end
        i=2;
        if t1==1&&t2==1&&t6==1
            alpha(:,2)=alpha(:,2)+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3))+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        end
        if t1==1&&t2==1&&t6==0
            alpha(:,2)=alpha(:,2)+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3));
        end
        if t1==1&&t2==0&&t6==1
            alpha(:,2)=alpha(:,2)+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        end
        if t1==0&&t2==1&&t6==1
            alpha(:,2)=alpha(:,2)+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3))+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        end
        if t1==1&&t2==0&&t6==0
            alpha(:,2)=alpha(:,2)+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1));
        end
        if t1==0&&t2==1&&t6==0
            alpha(:,2)=alpha(:,2)+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3));
        end
        if t1==0&&t2==0&&t6==1
            alpha(:,2)=alpha(:,2)+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        end
        if t1==0&&t2==0&&t6==0
            alpha(:,2)=alpha(:,2);
        end
        
        i=3;
        if t2==1&&t3==1
            alpha(:,3)=alpha(:,3)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        end
        if t2==1&&t3==0
            alpha(:,3)=alpha(:,3)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2));
        end
        if t2==0&&t3==1
            alpha(:,3)=alpha(:,3)+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        end
        if t2==0&&t3==0
            alpha(:,3)=alpha(:,3);
        end
        i=4;
        if t3==1&&t4==1&&t6==1
            alpha(:,4)=alpha(:,4)+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3));
        end
        if t3==1&&t4==1&&t6==0
            alpha(:,4)=alpha(:,4)+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3));
        end
        if t3==1&&t4==0&&t6==1
            alpha(:,4)=alpha(:,4)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3));
        end
        if t3==0&&t4==1&&t6==1
            alpha(:,4)=alpha(:,4)+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2));
        end
        if t3==1&&t4==0&&t6==0
            alpha(:,4)=alpha(:,4)+rho(i,k)*rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3));
        end
        if t3==0&&t4==1&&t6==0
            alpha(:,4)=alpha(:,4)+rho(i,k)*rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5));
        end
        if t3==0&&t4==0&&t6==1
            alpha(:,4)=alpha(:,4)+rho(i,k)*rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2));
        end
        if t3==0&&t4==0&&t6==0
            alpha(:,4)=alpha(:,4);
        end
        i=5;
        if t4==1&&t5==1
            alpha(:,5)=alpha(:,5)+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4))+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1));
        end
        if t4==1&&t5==0
            alpha(:,5)=alpha(:,5)+rho(i,k)*rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4));
        end
        if t4==0&&t5==1
            alpha(:,5)=alpha(:,5)+rho(i,k)*rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1));
        end
        if t4==0&&t5==0
            alpha(:,5)=alpha(:,5);
        end
        
        epsilon=0.99;
        i=1;
        if t1==1&&t5==1
            if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5)),inf)
                rho(1,k+1)=max([rho(1,k),rho(2,k),rho(5,k)]);
            else
                rho(1,k+1)=1.0001*max([rho(1,k),rho(2,k),rho(5,k)]);
            end
        end
        if t1==1&&t5==0
            if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2)),inf)
                rho(1,k+1)=max([rho(1,k),rho(2,k)]);
            else
                rho(1,k+1)=1.0001*max([rho(1,k),rho(2,k)]);
            end
        end
        if t1==0&&t5==1
            if norm(rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5)),inf)<=epsilon*norm(rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5)),inf)
                rho(1,k+1)=max([rho(1,k),rho(5,k)]);
            else
                rho(1,k+1)=1.0001*max([rho(1,k),rho(5,k)]);
            end
        end
        if t1==0&&t5==0
            rho(1,k+1)=rho(1,k);
        end
        
        
        i=2;
        if t1==1&&t2==1&&t6==1
            if norm(rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1))+rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3))+rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
                rho(2,k+1)=max([rho(1,k),rho(2,k),rho(3,k),rho(4,k)]);
            else
                rho(2,k+1)=1.0001*max([rho(1,k),rho(2,k),rho(3,k),rho(4,k)]);
            end
        end
        if t1==1&&t2==1&&t6==0
            if norm(rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3)),inf)<=epsilon*norm(rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1))+rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3)),inf)
                rho(2,k+1)=max([rho(1,k),rho(2,k),rho(3,k)]);
            else
                rho(2,k+1)=1.0001*max([rho(1,k),rho(2,k),rho(3,k)]);
            end
        end
        if t1==1&&t2==0&&t6==1
            if norm(rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1))+rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
                rho(2,k+1)=max([rho(1,k),rho(2,k),rho(4,k)]);
            else
                rho(2,k+1)=1.0001*max([rho(1,k),rho(2,k),rho(4,k)]);
            end
        end
        if t1==0&&t2==1&&t6==1
            if norm(rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3))+rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
                rho(2,k+1)=max([rho(2,k),rho(3,k),rho(4,k)]);
            else
                rho(2,k+1)=1.0001*max([rho(2,k),rho(3,k),rho(4,k)]);
            end
        end
        if t1==1&&t2==0&&t6==0
            if norm(rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1)),inf)<=epsilon*norm(rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1)),inf)
                rho(2,k+1)=max([rho(1,k),rho(2,k)]);
            else
                rho(2,k+1)=1.0001*max([rho(1,k),rho(2,k)]);
            end
        end
        if t1==0&&t2==1&&t6==0
            if norm(rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3)),inf)<=epsilon*norm(rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3)),inf)
                rho(2,k+1)=max([rho(2,k),rho(3,k)]);
            else
                rho(2,k+1)=1.0001*max([rho(2,k),rho(3,k)]);
            end
        end
        
        if t1==0&&t2==0&&t6==1
            if norm(rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
                rho(2,k+1)=max([rho(2,k),rho(4,k)]);
            else
                rho(2,k+1)=1.0001*max([rho(2,k),rho(4,k)]);
            end
        end
        if t1==0&&t2==0&&t6==0
            rho(2,k+1)=rho(2,k);
        end
        
        
        i=3;
        if t2==1&&t3==1
            if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
                rho(3,k+1)=max([rho(2,k),rho(3,k),rho(4,k)]);
            else
                rho(3,k+1)=1.0001*max([rho(2,k),rho(3,k),rho(4,k)]);
            end
        end
        if t2==1&&t3==0
            if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2)),inf)
                rho(3,k+1)=max([rho(2,k),rho(3,k)]);
            else
                rho(3,k+1)=1.0001*max([rho(2,k),rho(3,k)]);
            end
        end
        
        if t2==0&&t3==1
            if norm(rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
                rho(3,k+1)=max([rho(3,k),rho(4,k)]);
            else
                rho(3,k+1)=1.0001*max([rho(3,k),rho(4,k)]);
            end
        end
        if t2==0&&t3==0
            rho(3,k+1)=rho(3,k);
        end
        
        
        i=4;
        if t3==1&&t4==1&&t6==1
            if norm(rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3)),inf)<=epsilon*norm(rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5))+rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3)),inf)
                rho(4,k+1)=max([rho(2,k),rho(3,k),rho(4,k),rho(5,k)]);
            else
                rho(4,k+1)=1.0001*max([rho(2,k),rho(3,k),rho(4,k),rho(5,k)]);
            end
        end
        if t3==1&&t4==1&&t6==0
            if norm(rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3)),inf)<=epsilon*norm(rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5))+rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3)),inf)
                rho(4,k+1)=max([rho(3,k),rho(4,k),rho(5,k)]);
            else
                rho(4,k+1)=1.0001*max([rho(3,k),rho(4,k),rho(5,k)]);
            end
        end
        if t3==1&&t4==0&&t6==1
            if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2))+rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2))+rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3)),inf)
                rho(4,k+1)=max([rho(2,k),rho(3,k),rho(4,k)]);
            else
                rho(4,k+1)=1.0001*max([rho(2,k),rho(3,k),rho(4,k)]);
            end
        end
        if t3==0&&t4==1&&t6==1
            if norm(rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5))+rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2)),inf)<=epsilon*norm(rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5))+rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2)),inf)
                rho(4,k+1)=max([rho(2,k),rho(4,k),rho(5,k)]);
            else
                rho(4,k+1)=1.0001*max([rho(2,k),rho(4,k),rho(5,k)]);
            end
        end
        if t3==1&&t4==0&&t6==0
            if norm(rho(3,k)/(rho(i,k)+rho(3,k))*(x(:,i)-x(:,3)),inf)<=epsilon*norm(rho0(3)/(rho0(i)+rho0(3))*(x0(:,i)-x0(:,3)),inf)
                rho(4,k+1)=max([rho(3,k),rho(4,k)]);
            else
                rho(4,k+1)=1.0001*max([rho(3,k),rho(4,k)]);
            end
        end
        if t3==0&&t4==1&&t6==0
            if norm(rho(5,k)/(rho(i,k)+rho(5,k))*(x(:,i)-x(:,5)),inf)<=epsilon*norm(rho0(5)/(rho0(i)+rho0(5))*(x0(:,i)-x0(:,5)),inf)
                rho(4,k+1)=max([rho(4,k),rho(5,k)]);
            else
                rho(4,k+1)=1.0001*max([rho(4,k),rho(5,k)]);
            end
        end
        if t3==0&&t4==0&&t6==1
            if norm(rho(2,k)/(rho(i,k)+rho(2,k))*(x(:,i)-x(:,2)),inf)<=epsilon*norm(rho0(2)/(rho0(i)+rho0(2))*(x0(:,i)-x0(:,2)),inf)
                rho(4,k+1)=max([rho(2,k),rho(4,k)]);
            else
                rho(4,k+1)=1.0001*max([rho(2,k),rho(4,k)]);
            end
        end
        
        if t3==0&&t4==0&&t6==0
            rho(4,k+1)=rho(4,k);
        end
        
        i=5;
        if t1==1&&t5==1
            if norm(rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4))+rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1)),inf)<=epsilon*norm(rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4))+rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1)),inf)
                rho(5,k+1)=max([rho(1,k),rho(4,k),rho(5,k)]);
            else
                rho(5,k+1)=1.0001*max([rho(1,k),rho(4,k),rho(5,k)]);
            end
        end
        if t1==1&&t5==0
            if norm(rho(4,k)/(rho(i,k)+rho(4,k))*(x(:,i)-x(:,4)),inf)<=epsilon*norm(rho0(4)/(rho0(i)+rho0(4))*(x0(:,i)-x0(:,4)),inf)
                rho(5,k+1)=max([rho(4,k),rho(5,k)]);
            else
                rho(5,k+1)=1.0001*max([rho(4,k),rho(5,k)]);
            end
        end
        if t1==0&&t5==1
            if norm(rho(1,k)/(rho(i,k)+rho(1,k))*(x(:,i)-x(:,1)),inf)<=epsilon*norm(rho0(1)/(rho0(i)+rho0(1))*(x0(:,i)-x0(:,1)),inf)
                rho(5,k+1)=max([rho(1,k),rho(5,k)]);
            else
                rho(5,k+1)=1.0001*max([rho(1,k),rho(5,k)]);
            end
        end
        
        if t1==0&&t5==0
            rho(5,k+1)=rho(5,k);
        end
        
        
        x0=x;
        
        
        
        if k==100
            Xmat1_100=[Xmat1_100 x(:,1)];
            Xmat2_100=[Xmat2_100 x(:,2)];
            Xmat3_100=[Xmat3_100 x(:,3)];
            Xmat4_100=[Xmat4_100 x(:,4)];
            Xmat5_100=[Xmat5_100 x(:,5)];
        end
        if k==10
            Xmat1_10=[Xmat1_10 x(:,1)];
            Xmat2_10=[Xmat2_10 x(:,2)];
            Xmat3_10=[Xmat3_10 x(:,3)];
            Xmat4_10=[Xmat4_10 x(:,4)];
            Xmat5_10=[Xmat5_10 x(:,5)];
        end
        if k==500
            Xmat1_500=[Xmat1_500 x(:,1)];
            Xmat2_500=[Xmat2_500 x(:,2)];
            Xmat3_500=[Xmat3_500 x(:,3)];
            Xmat4_500=[Xmat4_500 x(:,4)];
            Xmat5_500=[Xmat5_500 x(:,5)];
        end
        if k==50
            Xmat1_50=[Xmat1_50 x(:,1)];
            Xmat2_50=[Xmat2_50 x(:,2)];
            Xmat3_50=[Xmat3_50 x(:,3)];
            Xmat4_50=[Xmat4_50 x(:,4)];
            Xmat5_50=[Xmat5_50 x(:,5)];
        end
        if k==1
            Xmat1_1=[Xmat1_1 x(:,1)];
            Xmat2_1=[Xmat2_1 x(:,2)];
            Xmat3_1=[Xmat3_1 x(:,3)];
            Xmat4_1=[Xmat4_1 x(:,4)];
            Xmat5_1=[Xmat5_1 x(:,5)];
        end
        if k==20
            Xmat1_20=[Xmat1_20 x(:,1)];
            Xmat2_20=[Xmat2_20 x(:,2)];
            Xmat3_20=[Xmat3_20 x(:,3)];
            Xmat4_20=[Xmat4_20 x(:,4)];
            Xmat5_20=[Xmat5_20 x(:,5)];
        end
        if k==30
            Xmat1_30=[Xmat1_30 x(:,1)];
            Xmat2_30=[Xmat2_30 x(:,2)];
            Xmat3_30=[Xmat3_30 x(:,3)];
            Xmat4_30=[Xmat4_30 x(:,4)];
            Xmat5_30=[Xmat5_30 x(:,5)];
        end
        if k==40
            Xmat1_40=[Xmat1_40 x(:,1)];
            Xmat2_40=[Xmat2_40 x(:,2)];
            Xmat3_40=[Xmat3_40 x(:,3)];
            Xmat4_40=[Xmat4_40 x(:,4)];
            Xmat5_40=[Xmat5_40 x(:,5)];
        end
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
