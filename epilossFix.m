
function NegativeValue=epilossFix(Para, Y_hat, Y_Obers)
ep=Para(1);
SD=Para(2);
res=Y_Obers-Y_hat';
l=length(Y_hat);
Llog=[];
for i=1:l
    if abs(res(i)/SD)<ep
        Llog(i)=0;
    else
       Llog(i)=abs(res(i)/SD)-ep;
    end
end
NegativeValue=sum(log(SD)+log(2*(1+ep))+Llog');
end