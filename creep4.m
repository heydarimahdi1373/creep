clc;
clear all;
close all;
T0=[3.42893	4.47099 4.48167	4.49771	4.50334	4.50895	4.51455	4.51999	4.52533	4.54103];
C0=[2.48E-11 3.00E-11	3.36E-11	3.92E-11	4.49E-11	5.02E-11	5.56E-11	5.87E-11	6.04E-11	6.15E-11];

TT=[0.001 0.008 0.064 0.512 4.096 32.768 262.164 2097.152 16777.216 134217.728]*10^-5;



T1=T0';
C1=C0';
n=numel(C0);
m=2;

%f1 = fit(T1(1:m),C1(1:m),'poly1');
%f2 = fit(T1(m:n),C1(m:n),'poly3');

%T(1:m)=(4.47099-3.42893)*rand(1,m)+3.42893*ones(1,m);
%clear rand
%T(m+1:n)=(4.54103-4.47099)*rand(1,n-m)+4.47099*ones(1,n-m);
%T=sort(T);
T=T0;
C=C0;

for i=1:m
    %C(i)=f1(T(i));
end
for i=m+1:n
    %C(i)=f2(T(i));
end

e=2.718281828;
X=-(e.^-(T'./TT))';
X(n,:)=1;
a=pinv(X')*C';

c=zeros(1,n);
for i=1:n
    for j=1:n-1
        c(i)=c(i)-a(j)*e^-(T(i)/TT(j));
    end
    c(i)=c(i)+a(n);
end

X=-(e.^-(T./TT'));
X(n,:)=1;
BB=zeros(n,1);
F=zeros(n,n);
for i=1:n
    for j=1:n
        for k=1:n
            F(i,j)=F(i,j)+X(i,k)*X(j,k);
        end
        BB(i)=BB(i)+X(i,j)*C(j);
    end
end
r=pinv(F)*BB;

c1=zeros(1,n);
for i=1:n
    for j=1:n-1
        c1(i)=c1(i)-r(j)*e^-(T(i)/TT(j));
    end
    c1(i)=c1(i)+r(n);
end

for i=1:n
    ax(i)=subplot(n,2,2*i-1);
    plot(ax(i),T,r(i)*e.^((T(i)*ones(1,n))./TT))
end
ax0=subplot(n,2,2*[1:n])
plot(T0,C0)
hold on
%plot(T(m:n),C(m:n))
%hold on
%plot(T,C)
%hold on
%plot(f1,T1,C1)
%hold on
%plot(f2,T0(m:n),C0(m:n))
%hold on
plot(T,c1)
hold on

