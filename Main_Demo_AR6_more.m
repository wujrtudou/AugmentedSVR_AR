%% This is the main code for our project "Augmented support vector regression with autoregressive process via iterative procedure".
%% Jinran Wu (email: jinran.wu@hdr.qut.edu.au) and You-Gan Wang (email: you-gam.wang@qut.edu.au)
%% School of Mathematical Sciences, Queensland Univerisity of Technology, Brisbane 4001, QLD, Australia
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TemporalSVR framework%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Clear the environment
%% Crude Oil price forecasting
%% Clear the environment
clc,clear
load crude_data
plot(data(1:716))
DataSet=[data(2:716,1),normalize([data(2:716,2:10),data(1:715,1)])];
%%
for i=1:1
InputTrain=DataSet(1:599+i,2:11);
OutputTrain=DataSet(1:599+i,1);
InputTest=DataSet(600+i,2:11);
OutputTest=DataSet(600+i,1);
%%
C=quantile(abs(OutputTrain),0.95);
g=0.01;
Epsilon=iqr(OutputTrain)/13.49;
[Alpha, Flag, B]=BasicSVR(InputTrain',OutputTrain', Epsilon, C,g);
BasicFitting=SVRPred(Alpha,Flag,B,InputTrain',g, InputTrain');
BasicPrediction(i)=SVRPred(Alpha,Flag,B,InputTrain',g, InputTest');
% Initialized series U
BasicU_Series=OutputTrain-BasicFitting';
% the temporal lag P
TemporalUInputTrain=[BasicU_Series(1:end-6)'; BasicU_Series(2:end-5)';BasicU_Series(3:end-4)';BasicU_Series(4:end-3)';BasicU_Series(5:end-2)';...
    BasicU_Series(6:end-1)'];
TemporalXInputTrain=InputTrain(7:end,:);
TemporalYOutTrain=OutputTrain(7:end);
% basic temporlSVR training
sigma=1;
[Alpha1, Flag1, B1]=TemporalSVR(TemporalXInputTrain',TemporalUInputTrain,TemporalYOutTrain',Epsilon,C,g,sigma);
% for iter=1:MaxIter
TemporalFitting1=TemporalSVRPred(Alpha1,Flag1,B1,TemporalXInputTrain',TemporalUInputTrain,g,TemporalXInputTrain',TemporalUInputTrain);
%figure(1)
%plot(TemporalXInputTrain,TrueTrain(2:end),'.k')
%hold on
%plot(TemporalXInputTrain,BasicFitting(2:end),'.r')
%hold on
%plot(TemporalXInputTrain,TemporalFitting1,'.b')
% Parameter estimation
x0=[1 1];
lb=[0 0];
ub=[];
A = [];
b = [];
Aeq = [];
beq = [];
ParameterEstor=fmincon(@(Para)epilossFix(Para, TemporalFitting1, TemporalYOutTrain),x0,A,b,Aeq,beq,lb,ub);
ep=ParameterEstor(1);
SD=ParameterEstor(2);
% Now we fix epsilon and scale to update omega, phi, and b.
SigmaValue=SD;
CT=quantile(abs(TemporalYOutTrain./SigmaValue),0.95);
EpsilonT=ep;
[Alpha2, Flag2, B2]=TemporalSVR(TemporalXInputTrain',TemporalUInputTrain,TemporalYOutTrain',EpsilonT,CT,g,SigmaValue);
% Now we update series U
BasicFitting3=SVRPred(Alpha2,Flag2,B2,TemporalXInputTrain',g, InputTrain');
NewU_Series=OutputTrain-BasicFitting3';
% Fix epsilon and scale, update omega, phi and b by new Useries
NewTemporalUInputTrain=[NewU_Series(1:end-6)';NewU_Series(2:end-5)';NewU_Series(3:end-4)';NewU_Series(4:end-3)';NewU_Series(5:end-2)'...
    ;NewU_Series(6:end-1)'];
[Alpha1, Flag1, B1]=TemporalSVR(TemporalXInputTrain',NewTemporalUInputTrain,TemporalYOutTrain',EpsilonT,CT,g,SigmaValue);
%% AIC calculation
K=6;
TemporalFittingFinal=TemporalSVRPred(Alpha1,Flag1,B1,TemporalXInputTrain',NewTemporalUInputTrain,g,TemporalXInputTrain',NewTemporalUInputTrain);
errors=TemporalYOutTrain-TemporalFittingFinal';
AIC_Value(i)=CalculateAIC(K,ep,SD,errors);
%end %  you can set the MaxIteration here
TemporalInputTest=NewU_Series(end);
TemporalPrediction(i)=TemporalSVRPred(Alpha1,Flag1,B1,TemporalXInputTrain',NewTemporalUInputTrain,g,InputTest',TemporalInputTest);
i
end
TrueResponse=DataSet(601:end,1)';
Temporal_MAE=mean(abs(TemporalPrediction-TrueResponse))
Basic_MAE=mean(abs(BasicPrediction-TrueResponse))
Temporal_RMAE=sqrt(mean((TemporalPrediction-TrueResponse).^2))
Basic_RMAE=sqrt(mean((BasicPrediction-TrueResponse).^2))
figure(1)
plot(TrueResponse)
hold on
plot(TemporalPrediction,'r')
hold on
plot(BasicPrediction,'g')
figure(2)
subplot(2,1,1)
parcorr(NewU_Series)
subplot(2,1,2)
autocorr(NewU_Series)