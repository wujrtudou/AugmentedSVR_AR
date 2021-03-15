function [Alpha, Flag, B]=TemporalSVR(X,U,Y,Epsilon,C,g,sigma)
l=size(X,2);
K1=[];
for i=1:l
    for j=1:l
        x1=X(:,i);
        y1=X(:,j);
        K1(i,j)=exp(-norm(x1-y1)^2*g);
    end
end
K2=[];
for i=1:l
    for j=1:l
        x2=U(:,i);
        y2=U(:,j);
        K2(i,j)=sum(x2.*y2);
    end
end

H=[K1, -K1; -K1, K1]+[K2, -K2; -K2, K2];

f=[Epsilon*sigma*ones(1,l)-Y, Epsilon*sigma*ones(1,l)+Y];
Aeq=[ones(1,l),-ones(1,l)];
Beq=0;
lb=0.*ones(2*l,1);
ub=C*sigma.*ones(2*l,1);
OPT=optimset;
OPT.LargeScale='off';
OPT.Display='off';
[Gamma,Obj]=quadprog(H,f,[],[],Aeq,Beq,lb,ub,[],OPT);
Alpha1=(Gamma(1:l,1))';
Alpha2=(Gamma((l+1):end,1))';
Alpha=Alpha1-Alpha2;
Flag=2*ones(1,l);
Err=0.0001;
for i=1:l
AA=Alpha1(i);
BB=Alpha2(i);
if (abs(AA-0)<=Err)&&(abs(BB-0)<=Err)
    Flag(i)=0;
end
if (AA>Err)&&(AA<C*sigma-Err)&&(abs(BB-0)<=Err)
Flag(i)=2;
end
if (abs(AA-0)<=Err)&&(BB>Err)&&(BB<C*sigma-Err)
Flag(i)=2;
end
if (abs(AA-C*sigma)<=Err)&&(abs(BB-0)<=Err)
Flag(i)=1;
end
if (abs(AA-0)<=Err)&&(abs(BB-C*sigma)<=Err)
Flag(i)=1;
end
end

B=0;
counter=0;
for i=1:l
AA=Alpha1(i);
BB=Alpha2(i);
if (AA>Err)&&(AA<C*sigma-Err)&&(abs(BB-0)<=Err)

SUM=0;
for j=1:l
if Flag(j)>0

SUM=SUM+Alpha(j)*(exp(-(norm(X(:,j)-X(:,i)))^2*g)+sum(U(:,j).*U(:,i)));
end
end
b=Y(i)-SUM-sigma*Epsilon;
B=B+b;
counter=counter+1;
end
if (abs(AA-0)<=Err)&&(BB>Err)&&(BB<C*sigma-Err)
SUM=0;
for j=1:l
if Flag(j)>0

SUM=SUM+Alpha(j)*(exp(-(norm(X(:,j)-X(:,i)))^2*g)+sum(U(:,j).*U(:,i)));
end
end
b=Y(i)-SUM+sigma*Epsilon;
B=B+b;
counter=counter+1;
end
end
B=B/counter;
end
