% Test Script for Iterative-Arnoldi-Tikhonov Method

% This test script requires that you have the IR-Tools toolbox downloaded
% and pointed to. You must also have the 'ArnoldiIterTik.m' code mapped to.

clear
clc

%% Problem details
N = 512; %PicSize
NoiseLevel = 0.01; 
options = PRset('trueImage', 'hst', 'BlurLevel', 'severe','BC','zero'); %reflexive,zero,periodic
[A,b,x,ProbInfo] = PRblurmotion(N,options);
[bn,NoiseInfo] = PRnoise(b,NoiseLevel);

%% IAT (iterative Arnoldi-Tikhonov Method)
x0 = zeros(N^2,1);
[Xapprox,RelResVec,RREvec,trackAlpha,trackExp,trackTikStep] = ArnoldiIterTik(A,bn,b,x,0.7,x0,NoiseLevel,80);

%% rrGMRES (via IRTools)

%Set-up options for rrGMRES soln
options = IRset('x_true',x,'NoiseLevel',NoiseLevel,'eta',1.01,'RegParam', 'discrep','MaxIter',300);
%use tau so that methods are comparable

%rrGMRES soln
[X,IterInfo] = IRrrgmres(A,bn,options);

%% Solution Info

method = ["rrGMRES";"IAT"];
finalRRE = [IterInfo.Enrm(length(IterInfo.Enrm));RREvec(length(RREvec))];
IterSteps = [length(IterInfo.Enrm);trackTikStep]; 
ArnoldiSteps = [length(IterInfo.Enrm);sum(trackExp)];

summaryTable = table(method,finalRRE,IterSteps,ArnoldiSteps)

%% Plotting 

subplot(2,3,1);
imagesc(reshape(bn,N,N))
axis image
axis off
colormap(gray)
title('Noised & Smoothed Data')

subplot(2,3,2);
imagesc(reshape(x,N,N))
axis image
axis off
colormap(gray)
title('True Image')

subplot(2,3,3);
imagesc(reshape(IterInfo.StopReg.X,N,N))
axis image
axis off
colormap(gray)
title('rrGMRES Soln (@breakout)')

subplot(2,3,6);
imagesc(reshape(Xapprox,N,N))
axis image
axis off
colormap(gray)
title('IAT Soln (@breakout)')

subplot(2,3,4)
breakout = (1.01*NoiseLevel).*ones(1,max(length(IterInfo.Rnrm),length(RelResVec)));
plot(1:length(breakout),breakout,'--k',1:length(IterInfo.Rnrm),IterInfo.Rnrm,'-og',1:length(RelResVec),RelResVec,'-*b')
title('Iteration vs. ||b^{\delta}-Ax_{approx}||/||b^{\delta}||')
xlabel('iteration #')
ylabel('||b^{\delta}-Ax_{approx}||/||b^{\delta}||')
legend('breakout_vec','rrGMRES','IAT')

subplot(2,3,5);
plot(1:length(IterInfo.Enrm),IterInfo.Enrm,'-og',1:length(RREvec),RREvec,'-*b')
hold on
title('Iteration vs. RRE')
legend('rrGMRES')
xlabel('iteration #')
ylabel('||x_{approx} - x||/||x||')
legend('rrGMRES','IAT')
