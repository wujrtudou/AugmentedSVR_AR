function Preds=TemporalSVRPred(Alpha,Flag,B,X,U,g,TestInput,TemporalInput)
l=length(Alpha);
s=size(TestInput,2);
Preds=[];
for j=1:s
OneInput1=TestInput(:,j);
OneInput2=TemporalInput(:,j);
SUM=0;
for i=1:l
if Flag(i)>0
SUM=SUM+Alpha(i)*(exp(-(norm(OneInput1-X(:,i)))^2*g)+sum(OneInput2.*U(:,i)));
end
Preds(j)=SUM+B;
end
end
end