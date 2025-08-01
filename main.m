%% clear
clc
clear
close all
%% Data Import
%OP is the input signal (Controller Output OP)
%fs is the sampling frequency (Hz)
%flag=1:oscillation;flag=0: No oscillation

load('benchmark.mat');
op=benchmark.chemicals.loop1.OP;
fs=1/benchmark.chemicals.loop1.Ts;
flag=benchmark.chemicals.loop1.oscillation;

% Ensure that the data is n * 1
[nx,ny]=size(op);
if nx<ny
    op=op';
end
%% Adaptive VMD
% The imf is a vector composed of kbest modalities for decomposition, denoted as N * kbest
% kbest is the optimal number of decompositions, indicating that the original signal has been decomposed into kbest modes
% kmax is the maximum number of modes 
[imf,kbest,kmax]=AdaptiveVMD(op,fs);
%% Select significant IMFs
% imf_s is significant IMFs
% sort is Modal original sorting
% n is number of imf_s
% delta is the growth rate of the correlation coefficient
[imf_s,sort,n,delta]=significant_IMFs(imf,op,kbest);
%% Oscillation probability assessment
%AD is a modal indicator value
AD=Osc(imf_s,n,fs,sort);