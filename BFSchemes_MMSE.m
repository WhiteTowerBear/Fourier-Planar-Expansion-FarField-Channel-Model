function [BFMatrix,MMSE]=BFSchemes_MMSE(ChannelPhase,var_noise)
 
%% %MMSE for Perfect channels
P =1/var_noise; %Linear scale
 

[U,S,V] = svd(ChannelPhase*ChannelPhase'+ var_noise * eye(size(ChannelPhase,1))); 
T=S;
T(find(S~=0)) = 1./S(find(S~=0));

BFMatrix =ChannelPhase'*  V * T' * U';
weights = ones(size(ChannelPhase,1),1);
 

BFMatrixDiagNorm=diag(sqrt(BFMatrix'*BFMatrix))';
BFMatrix = BFMatrix./repmat(BFMatrixDiagNorm,size(BFMatrix,1),1);
 
rhos = diag(abs(ChannelPhase*BFMatrix).^2)';
powerAllocationwMMSE = functionHeuristicPowerAllocation(rhos,P,weights);
 
BFMatrix = kron(sqrt(powerAllocationwMMSE),ones(size(BFMatrix,1),1)).*BFMatrix;


H_G_MMSE=ChannelPhase*BFMatrix;
MMSE_channelGains = abs(H_G_MMSE).^2;
MMSE_signalGains = diag(MMSE_channelGains);
MMSE_interferenceGains = sum(MMSE_channelGains,2)-MMSE_signalGains;
MMSE_rates = log2(1+MMSE_signalGains./(MMSE_interferenceGains+var_noise));
% MMSE_rates = log2(1+powerAllocationwMMSE./(var_noise));

MMSE=sum(MMSE_rates);
end