function Preds=SVRPred(Alpha,Flag,B,X,g,TestInput)
l=length(Alpha);
s=size(TestInput,2);
Preds=[];
for j=1:s
OneInput=TestInput(:,j);
SUM=0;
for i=1:l
if Flag(i)>0
SUM=SUM+Alpha(i)*exp(-(norm(OneInput-X(:,i)))^2*g);
end
Preds(j)=SUM+B;
end
end
end