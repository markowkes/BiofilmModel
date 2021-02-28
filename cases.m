function [inputvariables]=cases()
%% This function inputs 'TC' different test cases of input variables and produces matrix of test inputs which can be used to call on different segments for each case
% Each input variable will have its own array of values for each case (A,B,C,D,E,F based on chronological order in the array), and
% this function will produce an output 'inputvariables' matrix that can be
% organized into different cases and be called upon. 
mumax=[20 20 2 20 20 20];
Km   =[3 3 3 3 3 3000];
Yxs  =[0.5 0.5 0.5 0.5 0.5 0.5];
V    =[0.1 0.1 0.1 0.1 0.1 0.1];
Q    =[1 1 1 50 1 1];
A    =[1 1 1 1 1 1];
Sin  =[25 25 25 25 25 25];
So   =[25 25 25 25 25 25];
xo   =[10 10 10 10 10 10];
Daq  =[4.0E-5 4.0E-5 4.0E-5 4.0E-5 4.0E-5 4.0E-5];
De   =[1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E-5];
Xb   =[20000 20000 20000 20000 20000 20000];
Lfo  =[5.0E-6 3.0E-4 5.0E-6 5.0E-6 5.0E-6 5.0E-6];
LL   =[1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.0E-4 1.0E-4];
Kdet =[1900 1900 1900 1900 190000 1900];

inputvariables=[mumax; Km; Yxs; V; Q; A; Sin; So; xo; Daq; De; Xb; Lfo; LL; Kdet];





end