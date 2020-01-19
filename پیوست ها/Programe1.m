clc;
clear all;
close all;
T=[3.42893	3.5956	3.76747	3.93416	4.06439	4.17897	4.30924	4.42392	4.47099	4.48167	4.49771	4.50334	4.50895	4.51455	4.51999	4.52533	4.54103];
C=[2.48E-11	2.48E-11	2.49E-11	2.51E-11	2.54E-11	2.54E-11	2.62E-11	2.75E-11	3.00E-11	3.36E-11	3.92E-11	4.49E-11	5.02E-11	5.56E-11	5.87E-11	6.04E-11	6.15E-11];
n=numel(C);
e=2.718281828;
T=e.^(T);

X=(e.^-(T./T'));
X(n,:)=1;
B=zeros(n,1);
F=zeros(n,n);
for i=1:n
    for j=1:n
        for k=1:n
            F(i,j)=F(i,j)+X(i,k)*X(j,k);
        end
        B(i)=B(i)+X(i,j)*C(j);
    end
end
r=pinv(F)*B;

t=linspace(30,93,100);
m=numel(t);
c=zeros(1,m);
for i=1:m
    for j=1:n-1
        c(i)=c(i)+r(j)*e^-(t(i)/T(j));
    end
    c(i)=c(i)+r(n);
end

plot(T,C)
hold on
plot(t,c)