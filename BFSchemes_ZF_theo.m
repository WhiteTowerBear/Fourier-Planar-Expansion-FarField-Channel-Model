function [BFMatrix,ZF_theo]=BFSchemes_ZF_theo(RecVarianceVec,TraVarianceVec,var_noise,UsersNum,ChannelSigma_single)

%% The first method
n_s=length(TraVarianceVec);
n_r=length(RecVarianceVec);
 
 
[U,S,V] = svd(ChannelSigma_single*ChannelSigma_single'); 
T=S;
T(find(S~=0)) = 1./S(find(S~=0));

ChannelInv =ChannelSigma_single'* V * T' * U';

% GramMatrix=ChannelSigma_single*ChannelSigma_single';
% BFMatrix=ChannelSigma_single'*NeumannApprox(GramMatrix,3);

weights = ones(size(ChannelSigma_single,1),1);
 

for i=1:n_r
    ChannelColNorm(1,i)=norm(ChannelInv(:,i));
    BFMatrix(:,i)=sqrt(1/n_r)*ChannelInv(:,i)/ChannelColNorm(1,i);
end
P=1/var_noise;
rhos = diag(abs(ChannelSigma_single*BFMatrix).^2)';
p_u = functionHeuristicPowerAllocation(rhos,P,weights);


sigma_s=mean(TraVarianceVec);

  
nominator=(n_s-n_r+1)*sigma_s*RecVarianceVec;
denominator=var_noise*n_r;
ZF_temp = log2(1+p_u'.*nominator./denominator);

 

ZF_theo=sum(ZF_temp);

BFMatrix=[];

end