function [param]=cases(num)
%% This function inputs the 'num'th test case from an arbitrary number of testcases and produces a structure that produces the correct respective test parameters
% Each input variable will have its own array of values for each case
% (A,B,C,D,E,F based on chronological order in the array).
% This function will produce an output 'param' structure that can allow for
% each respective case to be called upon depending on the input num 

%Constants
mumax=[20 20 2 20 20 20 20 2000];
Km   =[3 3 3 3 3 3000 3 2500];
Yxs  =[0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
V    =[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
Q    =[1 1 1 50 1 1 1 1];
A    =[1 1 1 1 1 1 1 1];
Sin  =[25 25 25 25 25 25 25 25];
So   =[25 25 25 25 25 25 25 25];
xo   =[10 10 10 10 10 10 10 10];
Daq  =[4.0E-5 4.0E-5 4.0E-5 4.0E-5 4.0E-5 4.0E-5 4.0E2 4.00E-5];
De   =[1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E2 1.00E-5];
Xb   =[20000 20000 20000 20000 20000 20000 20000 20000];
Lfo  =[5.0E-6 3.0E-4 5.0E-6 5.0E-6 5.0E-6 5.0E-6 5.0E-6 5.00E-6];
LL   =[1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.00E-7];
Kdet =[1900 1900 1900 1900 190000 1900 1900 1900];

%Tank Parameters + Geometry
L=0.5; %[m]
W=0.5; %[m]
H=0.4; %[m]
SA=(V(num)/H)+2*((V(num)/L)+(V(num)/W)); %tank surface area [m^2] 

%Index variables under structure "param"
param.mumax=mumax(num);
param.Km   =Km(num);
param.Yxs  =Yxs(num);
param.V    =V(num);
param.Q    =Q(num);
param.A    =A(num);
param.Sin  =Sin(num);
param.So   =So(num);
param.xo   =xo(num);
param.Daq  =Daq(num);
param.De   =De(num);
param.Xb   =Xb(num);
param.Lfo  =Lfo(num);
param.LL   =LL(num);
param.Kdet =Kdet(num);
param.SA   =SA;
end