function [BFMatrix,ZF]=BFSchemes_ZF_noinverse(ChannelPhase,var_noise,ChannelRandom,RecVarianceVec,TraVarianceVec)

%% ZF for Perfect channels
P =1/var_noise; %Linear scale 
% GramMatrix=ChannelRandom*ChannelRandom';
% WW_inverse=NeumannApprox(GramMatrix,5);
% W_inverse=ChannelRandom'*WW_inverse;
% HaHa_inverse=diag(1./RecVarianceVec)*W_inverse'*diag(1./TraVarianceVec.^2) ...
%     *W_inverse*(diag(1./RecVarianceVec));

GramMatrix=ChannelRandom*diag(TraVarianceVec.^2)*ChannelRandom';
GramMatrixNorm=norm(GramMatrix,'fro');
GramMatrix=GramMatrix/GramMatrixNorm;
WW_inverse=NeumannApprox(GramMatrix,4)/GramMatrixNorm;
HaHa_inverse=diag(1./RecVarianceVec)*WW_inverse*diag(1./RecVarianceVec);
ChannelInv =ChannelPhase'* HaHa_inverse; 
 

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