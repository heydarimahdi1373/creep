clc;
clear all;
close all;
T0=[3.42893	3.5956	3.76747	3.93416	4.06439	4.17897	4.30924	4.42392	4.47099	4.48167	4.49771	4.50334	4.50895	4.51455	4.51999	4.52533	4.54103];
C0=[2.48E-11	2.48E-11	2.49E-11	2.51E-11	2.54E-11	2.54E-11	2.62E-11	2.75E-11	3.00E-11	3.36E-11	3.92E-11	4.49E-11	5.02E-11	5.56E-11	5.87E-11	6.04E-11	6.15E-11];
T1=T0';
C1=C0';
n=numel(C0);
m=9;

%T1=T0';
%C1=C0';
f1 = fit(T1(1:m),C1(1:m),'poly4');
%T2=T1(m:n);
%C2=C1(m:n);
f2 = fit(T1(m:n),C1(m:n),'poly3');

%p1=[4.765E-11 -7.357E-10 4.252E-09 -1.09E-08 1.048E-08];
%q1=@(x) p1(1)*x^4+p1(2)*x^3+p1(3)*x^2+p1(4)*x+p1(5);
%p2=[-0.008163 0.2212 -2.498 15.04 -50.95 92.04 -69.28];
%q2=@(x) p2(1)*x^6+p2(2)*x^5+p2(3)*x^4+p2(4)*x^3+p2(5)*x^2+p2(6)*x+p2(7);

T(1:m)=(4.47099-3.42893)*rand(1,m)+3.42893*ones(1,m);
clear rand
T(m+1:n)=(4.54103-4.47099)*rand(1,n-m)+4.47099*ones(1,n-m);
T=sort(T);


for i=1:m
    %C(i)=q1(T(i));
    C(i)=f1(T(i));
end
for i=m+1:n
    %C(i)=q2(T(i));
    C(i)=f2(T(i));
end

e=2.718281828;
X=(e.^-(T'./T))';
X(n,:)=1;
a=pinv(X')*C';

c=zeros(1,n);
for i=1:n
    for j=1:n-1
        c(i)=c(i)+a(j)*e^-(T(i)/T(j));
    end
    c(i)=c(i)+a(n);
end

X=(e.^-(T./T'));
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
        c1(i)=c1(i)+r(j)*e^-(T(i)/T(j));
    end
    c1(i)=c1(i)+r(17);
end

%T1=T0';
%C1=C0';
%f1 = fit(T1(1:m),C1(1:m),'poly4');
%T2=T1(m:n);
%C2=C1(m:n);
%f2 = fit(T2,C2,'poly3');
TT=T';

plot(T0,C0)
hold on
%plot(T(m:n),C(m:n))
hold on
plot(T,c1,T,c)
hold on
%plot(f1,T1,C1)
%hold on
%plot(f2,T0(m:n),C0(m:n))
%hold on