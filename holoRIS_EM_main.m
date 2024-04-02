clc
clear all

global FRAME
FRAME=200;
double=0; % the double scattering environment 
step=2;
SNR_Range=-10:step:20;

% global spacing % the spacing between RIS elements
global lambda
global RecLength
global TraLength
lambda=0.03; % wavelength 
UsersNum=1; 
UserSpacing=100*lambda;
RecSpacing_num=3;
RecSpacing=lambda/RecSpacing_num;
TraSpacing_num=3;
TraSpacing=lambda/TraSpacing_num;

Nr_X=12; Nr_Y=12;
Ns_X=27; Ns_Y=27;

RecLength.L_x=Nr_X*RecSpacing;RecLength.L_y=Nr_Y*RecSpacing;
RecLength.distance=[6,4,2];

TraLength.L_x=Ns_X*TraSpacing;TraLength.L_y=Ns_Y*TraSpacing;
TraLength.distance=[0,0,0];
% % The length of received  surface
% RecLength =struct('L_x',2*lambda,'L_y',4*lambda,'distance',[6,4,2]);
% % The length of transmitted  surface
% TraLength=struct('L_x',16*lambda,'L_y',10*lambda,'distance',[0,0,0]);

%% Generate the channel matrix
% Obtain the sampled eclipse and their coordinates given a surface
[RecSamplePoint,WD_RecX_vec,WD_RecY_vec]=WD_eclipse_sample(lambda,RecLength);
[TraSamplePoint,WD_TraX_vec,WD_TraY_vec]=WD_eclipse_sample(lambda,TraLength);

[RecVarianceVec,WD_RecX_vec,WD_RecY_vec,nr_act]=WaveDomainChannel(RecSamplePoint,WD_RecX_vec,WD_RecY_vec,lambda,RecLength);
[TraVarianceVec,WD_TraX_vec,WD_TraY_vec,ns_act]=WaveDomainChannel(TraSamplePoint,WD_TraX_vec,WD_TraY_vec,lambda,TraLength);

% average
if double~=1
    tt1=mean(TraVarianceVec);
    TraVarianceVec=tt1*ones(ns_act,1);
end
ChannelSigma_single=RecVarianceVec*TraVarianceVec';
if UsersNum>1
    ChannelSigma=repmat(ChannelSigma_single,UsersNum,1);
    RecVarianceVec=repmat(RecVarianceVec,UsersNum,1);
else
    ChannelSigma=ChannelSigma_single;
end

[TraResVector,TraNumNs]=ResponseGenerate(TraSpacing,lambda,TraLength,WD_TraX_vec,WD_TraY_vec,'Transmit'); 


% Multi-user Case
RecResVector=[];
for user=1:UsersNum 
    RecLength(user)=RecLength(1);
    RecLength(user).distance(1) =(user-1)*UserSpacing+RecLength(1).distance(1);
    tempUserInfo=RecLength(user);
    [RecResVectorPart,RecNumNr]=ResponseGenerate(RecSpacing,lambda,tempUserInfo,WD_RecX_vec,WD_RecY_vec,'Receive');
    RecResVector=blkdiag(RecResVector,RecResVectorPart);
end

%% Pass Channel 
% [Phi]=PhaseGenerate(TraNumNs);  % Phase matrix 
Phi=eye(TraNumNs); % Phi is an identity matrix
TransSig=randn(UsersNum*RecNumNr,1)+1i*randn(UsersNum*RecNumNr,1); % transmitted signal

MRT_snr=[];
ZF_snr=[];
ZF_NS_snr_theo=[];
MMSE_snr=[];
MRT_snr_theo=[];
ZF_snr_theo=[];
MMSE_snr_theo=[];
for SNR=SNR_Range(1):step:SNR_Range(end)
    var_noise=10^(-0.1*SNR);
    MRT_collect=[];
    ZF_collect=[];
    MMSE_collect=[];
    ZF_NS_collect=[];
    
    MRT_theo_collect=[];
    ZF_iter_collect=[];
    MMSE_theo_collect=[];
    alpha_mrt_collect=[];
 
    for frame=1:FRAME
        noise=sqrt(var_noise/2)*(randn(UsersNum*RecNumNr,1)+1i*randn(UsersNum*RecNumNr,1));
        ChannelRandom=sqrt(1/2)*(randn(UsersNum*nr_act,ns_act)+1i*randn(UsersNum*nr_act,ns_act));
              
%         ChannelHa=ChannelSigma.*ChannelRandom;
        ChannelHa= diag(RecVarianceVec)*ChannelRandom*diag(TraVarianceVec);
%         Channel=1/sqrt(UsersNum*RecNumNr*TraNumNs)*RecResVector*ChannelHa*TraResVector;
        Channel= ChannelHa*TraResVector;
        ChannelPhase=Channel*Phi;
        % System Model
        [BFMatrix_MMSE,MMSE_iter]=BFSchemes_MMSE (ChannelHa,var_noise);
        MMSE_collect=[MMSE_collect,MMSE_iter];       
        [BFMatrix_ZF,ZF_iter]=BFSchemes_ZF (ChannelHa,var_noise);
        ZF_collect=[ZF_collect,ZF_iter];      
        [~,ZF_NS_iter]=BFSchemes_ZF_noinverse(ChannelHa,var_noise,ChannelRandom,RecVarianceVec,TraVarianceVec);
        ZF_NS_collect=[ZF_NS_collect,ZF_NS_iter]; 
        [BFMatrix_MRT,MRT_iter]=BFSchemes_MRT(ChannelHa,var_noise);
        MRT_collect=[MRT_collect,MRT_iter];      
        
    end
    
%     [BFMatrixMMSE_theo,MMSE_iter_theo]=BFSchemes_MMSE_theo(RecVarianceVec,TraVarianceVec,var_noise,UsersNum,ChannelSigma);
    [BFMatrix_ZF_theo,ZF_iter_theo]=BFSchemes_ZF_theo(RecVarianceVec,TraVarianceVec,var_noise,UsersNum,ChannelSigma);
    [BFMatrixMRT_theo,MRT_iter_theo]=BFSchemes_MRT_theo(RecVarianceVec,TraVarianceVec,var_noise,UsersNum,ChannelSigma);
 
    MRT_snr=[MRT_snr,sum(MRT_collect)/FRAME];
    ZF_snr=[ZF_snr,sum(ZF_collect)/FRAME];
    ZF_NS_snr_theo=[ZF_NS_snr_theo,sum(ZF_NS_collect)/FRAME ];
    MMSE_snr=[MMSE_snr,sum(MMSE_collect)/FRAME];
   
%     MMSE_snr_theo=[MMSE_snr_theo,MMSE_iter_theo];
    ZF_snr_theo=[ZF_snr_theo,ZF_iter_theo ];
    MRT_snr_theo=[MRT_snr_theo,MRT_iter_theo ]; 
    
    
    
end

figure
plot(SNR_Range,MRT_snr,'s-m','linewidth',2)
hold on
plot(SNR_Range,ZF_snr,'o-b','linewidth',2)
hold on 
plot(SNR_Range,ZF_NS_snr_theo,'x-k','linewidth',2)
hold on 
plot(SNR_Range,MMSE_snr,'d-r','linewidth',2)
hold on
plot(SNR_Range,MRT_snr_theo,'s:m','linewidth',2)
hold on
plot(SNR_Range,ZF_snr_theo,'o:b','linewidth',2)
hold on 
% plot(SNR_Range,MMSE_snr_theo,'d:r','linewidth',2)
% legend('MRT ','ZF ','MMSE')
% legend('MRT (theo)','ZF (theo)','MMSE (theo)')
legend('MRT ','ZF ','ZF (NS)','MMSE' ,'MRT (theo)','ZF (theo)' ,'Interpreter','latex' )

xlabel('SNR (dB)','Interpreter','latex')
ylabel('Spectral Efficiency (bits/s/Hz)','Interpreter','latex')
grid on
title(['$N_s$=',num2str(TraNumNs),', $N_r$=',num2str(RecNumNr), ...
    ', $\Delta_s$=$\lambda$/',num2str(TraSpacing_num), ...
    ', $\Delta_r$=$\lambda$/',num2str(RecSpacing_num)],'Interpreter','latex')

  
