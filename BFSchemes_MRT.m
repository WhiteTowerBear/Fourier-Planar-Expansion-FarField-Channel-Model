function [BFMatrix,MRT]=BFSchemes_MRT(ChannelPhase,var_noise)

%% MRC for Perfect channels
P =1/var_noise; %Linear scale


BFMatrix=ChannelPhase';
weights = ones(size(BFMatrix,2),1);
 
 
BFMatrix=BFMatrix/norm(BFMatrix,'fro');

H_G_MRC=ChannelPhase*BFMatrix;
MRC_signalGains = diag(abs(H_G_MRC).^2)';
powerAllocationMRT = functionHeuristicPowerAllocation(MRC_signalGains,P,weights);
BFMatrix = kron(sqrt(powerAllocationMRT),ones(size(BFMatrix,1),1)).*BFMatrix;

channelGains = abs(ChannelPhase*BFMatrix).^2;
signalGains = diag(channelGains);
interferenceGains = sum(channelGains,2)-signalGains;
rates = log2(1+signalGains./(interferenceGains+var_noise  ));

 MRT = sum(rates);

 

end