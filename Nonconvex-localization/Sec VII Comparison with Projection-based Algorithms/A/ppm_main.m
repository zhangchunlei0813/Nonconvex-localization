clc,clear all,close all,
load 'rmat_nf001.mat'
load 'X_initial_nf001.mat'
m=5; %anchors number
D=2; %D-dimension
L=1000;

a=[5.2 8.9 12.1 10.5 6.3;1.7 5.0 9.0 13.7 12.2];
target=[7.5;7.5];

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

XX=zeros(2,m);
for l=1:L
    r=rmat(1+5*(l-1):5+5*(l-1));
    x0(:,1)=X_initial(:,l);
    tic
    for j=1:kk
        x=0;
        for i=1:m
            x=x+a(:,i)+ r(i)*(x0(:,j)-a(:,i))/norm(x0(:,j)-a(:,i));
        end
        x=x/m;
        for i=1:m
            XX(:,i)=a(:,i)+ r(i)*(x0(:,j)-a(:,i))/norm(x0(:,j)-a(:,i));
        end
        x0(:,j+1)=x;
        if j==101
            Xmat100=[Xmat100 x];
        end
        if j==11
            Xmat10=[Xmat10 x];
        end
        if j==501
            Xmat500=[Xmat500 x];
        end
        if j==51
            Xmat50=[Xmat50 x];
        end
        if j==2
            Xmat1=[Xmat1 x];
        end
        if j==21
            Xmat20=[Xmat20 x];
        end
        if j==31
            Xmat30=[Xmat30 x];
        end
        if j==41
            Xmat40=[Xmat40 x];
        end
    end
    toc
    Xmat=[Xmat x];
end

Xmat_sum=0;
for nn=1:L
    Xmat_sum=Xmat_sum+(norm(Xmat(:,nn)-target))^2;
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



