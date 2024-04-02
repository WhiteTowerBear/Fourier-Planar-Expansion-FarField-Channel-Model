% This is for SPM simulation
clc
clear all

global FRAME
FRAME=200;
double=1; % the double scattering environment 
 

% global spacing % the spacing between RIS elements
global lambda
global RecLength
global TraLength
lambda=1; % wavelength 
UsersNum=1; 
UserSpacing=100*lambda;
RecSpacing_num=3;
RecSpacing=lambda/RecSpacing_num;
TraSpacing_num=15;
TraSpacing=lambda/TraSpacing_num;

Nr_X=24; Nr_Y=24;
Nr=Nr_X*Nr_Y;
Ns_X=30; Ns_Y=30;
Ns=Ns_X*Ns_Y;


RecLength.L_x=Nr_X*RecSpacing;RecLength.L_y=Nr_Y*RecSpacing;
RecLength.distance=[0,0,10];

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
 
R_lambda_r=diag(reshape(RecVarianceVec.*RecVarianceVec,[],1));
 R_Cor_r=RecResVector*R_lambda_r*RecResVector';


 
eigenvalues_r= sort(abs(eig(R_Cor_r)),'descend');

% R_lambda_s=diag(reshape(TraVarianceVec.*TraVarianceVec,[],1));
%  R_Cor_s=TraResVector'*R_lambda_s*TraResVector;
% eigenvalues_s= sort(eig(R_Cor_s),'descend');

figure;
hold on; box on;
plot(1:Nr,eigenvalues_r,'b-','LineWidth',2);
% plot(1:Ns,eigenvalues_s,'k--','LineWidth',2);
set(gca,'Yscale','log');
ylim([1e-2 1e2]);
legend({'$Rec$','$Tran$'} ,'Interpreter','latex','Location','NorthEast');


xlabel('Eigenvalue number','Interpreter','latex');
ylabel('Eigenvalue','Interpreter','latex');
grid on
title(['$N_s$=',num2str(TraNumNs),', $N_r$=',num2str(RecNumNr), ...
    ', $\Delta_s$=$\lambda$/',num2str(TraSpacing_num), ...
    ', $\Delta_r$=$\lambda$/',num2str(RecSpacing_num)],'Interpreter','latex')

 

 
 % Compute channel capacity
 iter=1;
 noise=1; 
snr_range=-5:2:15;
 for SNR_value=snr_range 
     trans_power=noise*10^(0.1*SNR_value);
     R_lambda_r(R_lambda_r==0)=[];
     Capacity(1,iter)=sum(log2(1+trans_power*R_lambda_r));
     
     iter=iter+1;
 end
 figure;
hold on;   box on;
plot(snr_range,Capacity,'ro-','LineWidth',2);
legend('Capacity' ,'Interpreter','latex','Location','NorthEast');
grid on
title(['$N_s$=',num2str(TraNumNs),', $N_r$=',num2str(RecNumNr), ...
    ', $\Delta_s$=$\lambda$/',num2str(TraSpacing_num), ...
    ', $\Delta_r$=$\lambda$/',num2str(RecSpacing_num)],'Interpreter','latex')

