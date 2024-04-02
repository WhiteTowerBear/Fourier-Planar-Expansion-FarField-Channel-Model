function [Phi]=PhaseGenerate(N)
% This function aims at generating random pilot phase matrix
 
Phi=diag(exp(1i*2*pi*rand(N,1)));% randomly generation

% plot(Phi,'.')
 
 