clc,clear all,close all,
load 'rmat_nf001.mat'
load 'X_initial_nf001.mat'
load 'NN_05.mat'
load 'TT_05.mat'

m=5; %anchors number
D=2; %D-dimension
L=1000;

a=[5.2 8.9 12.1 10.5 6.3;1.7 5.0 9.0 13.7 12.2];
target=[7.5;7.5];
d=zeros(1,m);
for i=1:m
    d(i)=norm(a(:,i)-target);
end

kk=1000;
Xmat=[];
Xmat1=[];
Xmat10=[];
Xmat20=[];
Xmat30=[];
Xmat40=[];
Xmat50=[];
Xmat100=[];
Xmat500=[];

BETA=[];
for j=1:L
    r=rmat(1+5*(j-1):5+5*(j-1));

    alpha=zeros(2,kk);

    alpha(:,1)=X_initial(:,j);
    for k=1:kk
        beta=alpha(:,k);
        if TT(k,1+6*(j-1):6*j-1)==[1 1 1 1 1]; 
            for i=1:m
                alpha(:,k+1)=a(:,i)+r(i)*(beta-a(:,i))/norm(beta-a(:,i));
                beta=alpha(:,k+1);
            end
        else
            alpha(:,k+1)=alpha(:,k);
            beta=alpha(:,k+1);
        end
        if k==100
            Xmat100=[Xmat100 beta];
        end
        if k==10
            Xmat10=[Xmat10 beta];
        end
        if k==500
            Xmat500=[Xmat500 beta];
        end
        if k==50
            Xmat50=[Xmat50 beta];
        end
        if k==1
            Xmat1=[Xmat1 beta];
        end
        if k==20
            Xmat20=[Xmat20 beta];
        end
        if k==30
            Xmat30=[Xmat30 beta];
        end
        if k==40
            Xmat40=[Xmat40 beta];
        end
    end
    BETA=[BETA beta];
end

Xmat_sum=0;
for nn=1:L
    Xmat_sum=Xmat_sum+(norm(BETA(:,nn)-target))^2;
end
Xmat_average=sqrt(Xmat_sum/L);

Xmat_sum10=0;
for nn=1:L
    Xmat_sum10=Xmat_sum10+(norm(Xmat10(:,nn)-target))^2;
end
Xmat_average10=sqrt(Xmat_sum10/L);

Xmat_sum1=0;
for nn=1:L
    Xmat_sum1=Xmat_sum1+(norm(Xmat1(:,nn)-target))^2;
end
Xmat_average1=sqrt(Xmat_sum1/L);

Xmat_sum100=0;
for nn=1:L
    Xmat_sum100=Xmat_sum100+(norm(Xmat100(:,nn)-target))^2;
end
Xmat_average100=sqrt(Xmat_sum100/L);

Xmat_sum500=0;
for nn=1:L
    Xmat_sum500=Xmat_sum500+(norm(Xmat500(:,nn)-target))^2;
end
Xmat_average500=sqrt(Xmat_sum500/L);

Xmat_sum50=0;
for nn=1:L
    Xmat_sum50=Xmat_sum50+(norm(Xmat50(:,nn)-target))^2;
end
Xmat_average50=sqrt(Xmat_sum50/L);

Xmat_sum20=0;
for nn=1:L
    Xmat_sum20=Xmat_sum20+(norm(Xmat20(:,nn)-target))^2;
end
Xmat_average20=sqrt(Xmat_sum20/L);

Xmat_sum30=0;
for nn=1:L
    Xmat_sum30=Xmat_sum30+(norm(Xmat30(:,nn)-target))^2;
end
Xmat_average30=sqrt(Xmat_sum30/L);

Xmat_sum40=0;
for nn=1:L
    Xmat_sum40=Xmat_sum40+(norm(Xmat40(:,nn)-target))^2;
end
Xmat_average40=sqrt(Xmat_sum40/L);


