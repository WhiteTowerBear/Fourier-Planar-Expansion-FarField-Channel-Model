function [BFMatrix,ZF]=BFSchemes_ZF(ChannelPhase,var_noise)

%% ZF for Perfect channels
P =1/var_noise; %Linear scale

[U,S,V] = svd(ChannelPhase*ChannelPhase'); 
T=S;
T(find(S~=0)) = 1./S(find(S~=0));
 
ChannelInv =ChannelPhase'* V * T' * U'; 
 

weights = ones(size(ChannelPhase,1),1);
[n_r,~]=size(ChannelPhase);

for i=1:n_r
    ChannelColNorm(1,i)=norm(ChannelInv(:,i));
    BFMatrix(:,i)=sqrt(1/n_r)*ChannelInv(:,i)/ChannelColNorm(1,i);
end

% BFMatrixDiagNorm=diag(sqrt(ChannelInv'*ChannelInv))';
% BFMatrix = sqrt(1/n_r)*ChannelInv./repmat(BFMatrixDiagNorm,size(ChannelInv,1),1);
 
rhos = diag(abs(ChannelPhase*BFMatrix).^2)';
powerAllocationwZFBF = functionHeuristicPowerAllocation(rhos,P,weights);
 
% ZF_rate =  log2(1+powerAllocationwZFBF./( var_noise ));

ZF_rate =  log2(1+powerAllocationwZFBF./(n_r*var_noise*ChannelColNorm.^2));
ZF=sum(ZF_rate);

end