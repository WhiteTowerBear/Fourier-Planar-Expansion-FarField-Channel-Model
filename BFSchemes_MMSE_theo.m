function [BFMatrix,MMSE_theo]=BFSchemes_MMSE_theo(RecVarianceVec,TraVarianceVec,var_noise,UsersNum,ChannelSigma_single)
 
%% %MMSE for Perfect channels
P =1/var_noise; %Linear scale
[U,S,V] = svd(ChannelSigma_single*ChannelSigma_single'+ var_noise * eye(size(ChannelSigma_single,1))); 
T=S;
T(find(S~=0)) = 1./S(find(S~=0));

BFMatrix =ChannelSigma_single'*  V * T' * U';
weights = ones(size(ChannelSigma_single,1),1);
BFMatrixDiagNorm=diag(sqrt(BFMatrix'*BFMatrix))';
BFMatrix = BFMatrix./repmat(BFMatrixDiagNorm,size(BFMatrix,1),1);
rhos = diag(abs(ChannelSigma_single*BFMatrix).^2)';
p_u = functionHeuristicPowerAllocation(rhos,P,weights);
 

n_s=length(TraVarianceVec);
n_r=length(RecVarianceVec);
sigma_s=mean(TraVarianceVec);  
nominator=var_noise^2./RecVarianceVec.^2./p_u';
denominator=p_u'*n_s*sigma_s^2.*RecVarianceVec.^2./(sum(1./RecVarianceVec.^2) ...
    - 1./RecVarianceVec.^2)+var_noise*n_r/(n_s-n_r)*(n_s^2+n_s)*sigma_s^2;

MMSE_temp = log2(1+ nominator./denominator);
MMSE_theo=sum(MMSE_temp);

BFMatrix=[];
end