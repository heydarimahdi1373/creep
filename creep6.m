clc;
clear all;
close all;
T=[30.84362241	36.43755596	43.27045199	51.1191918	58.22937777	65.29856105	74.38393566	83.42266208	87.44324899	88.38214769	89.81122792	90.31829118	90.8264007	91.33645537	91.83467963	92.32638851	93.78735136];
C=[2.48E-11	2.48E-11	2.49E-11	2.51E-11	2.54E-11	2.54E-11	2.62E-11	2.75E-11	3.00E-11	3.36E-11	3.92E-11	4.49E-11	5.02E-11	5.56E-11	5.87E-11	6.04E-11	6.15E-11]*10^11;
e=2.718281828;
n=numel(C);
t=linspace(0,400,100);
tt=1E2;
mm=5/(14*tt);
k=6.236E-12/mm;
g=1-e.^(-t/tt)

f=100*(1-e.^(-(t-85.5)/tt));
g=e.^(-(t)/tt);
plot(T,C,t,f)
plot(t,g)