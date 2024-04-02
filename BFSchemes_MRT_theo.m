function [BFMatrix,MRT]=BFSchemes_MRT_theo(RecVarianceVec,TraVarianceVec,var_noise,UsersNum,ChannelSigma,alpha_mrt)

%% MRC for Perfect channels
P =1/var_noise; %Linear scale
BFMatrix=ChannelSigma';
weights = ones(size(BFMatrix,2),1);
BFMatrix=BFMatrix/norm(BFMatrix,'fro');
H_G_MRC=ChannelSigma*BFMatrix;
MRC_signalGains = diag(abs(H_G_MRC).^2)';
p_u = functionHeuristicPowerAllocation(MRC_signalGains,P,weights);

n_s=length(TraVarianceVec);
n_r=length(RecVarianceVec);
sigma_s=mean(TraVarianceVec);
TraVarianceVec=TraVarianceVec.^2;
RecVarianceVec=RecVarianceVec.^2;
  
% nominator=p_u*(n_s-1)*(n_s-2) * sigma_s * RecVarianceVec.^2;
% denominator=p_u*(n_s-2) *sigma_s*RecVarianceVec.*(sum(RecVarianceVec)-RecVarianceVec) ...
%     +var_noise*n_s*sum(RecVarianceVec);

% For the case M>>2
nominator=p_u'*n_s * sigma_s .* RecVarianceVec.^2;
denominator=p_u'*sigma_s.*RecVarianceVec.*(sum(RecVarianceVec)-RecVarianceVec) ...
    +var_noise*sum(RecVarianceVec);
  
%%
MRT_temp = log2(1+nominator./denominator);
MRT=sum(MRT_temp);

BFMatrix=[];
end