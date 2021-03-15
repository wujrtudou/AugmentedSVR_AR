%% AIC calculatin
function AIC_Value=CalculateAIC(K,eps,scale,error)
for i=1:length(error)
    k=error(i)/scale;
   if k>eps
       l1(i)=k-eps;
   elseif k<-eps
           l1(i)=-k-eps;
   else
       l1(i)=0;
   end
end
loss=sum(l1);
a1=2*K;
a2=length(error)*log(scale)+length(error)*log(2*(1+eps))+loss;
AIC_Value=a1+2*a2;
end