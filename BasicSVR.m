function [Alpha, Flag, B]=BasicSVR(X,Y,Epsilon,C,g)
l=size(X,2);
K=[];
for i=1:l
    for j=1:l
        x=X(:,i);
        y=X(:,j);
        K(i,j)=exp(-norm(x-y)^2*g);
    end
end
H=[K, -K; -K, K];
f=[Epsilon*ones(1,l)-Y, Epsilon*ones(1,l)+Y];
Aeq=[ones(1,l),-ones(1,l)];
Beq=0;
lb=0.*ones(2*l,1);
ub=C.*ones(2*l,1);
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
if (AA>Err)&&(AA<C-Err)&&(abs(BB-0)<=Err)
Flag(i)=2;
end
if (abs(AA-0)<=Err)&&(BB>Err)&&(BB<C-Err)
Flag(i)=2;
end
if (abs(AA-C)<=Err)&&(abs(BB-0)<=Err)
Flag(i)=1;
end
if (abs(AA-0)<=Err)&&(abs(BB-C)<=Err)
Flag(i)=1;
end
end

B=0;
counter=0;
for i=1:l
AA=Alpha1(i);
BB=Alpha2(i);
if (AA>Err)&&(AA<C-Err)&&(abs(BB-0)<=Err)

SUM=0;
for j=1:l
if Flag(j)>0

SUM=SUM+Alpha(j)*exp(-(norm(X(:,j)-X(:,i)))^2*g);
end
end
b=Y(i)-SUM-Epsilon;
B=B+b;
counter=counter+1;
end
if (abs(AA-0)<=Err)&&(BB>Err)&&(BB<C-Err)
SUM=0;
for j=1:l
if Flag(j)>0

SUM=SUM+Alpha(j)*exp(-(norm(X(:,j)-X(:,i)))^2*g);
end
end
b=Y(i)-SUM+Epsilon;
B=B+b;
counter=counter+1;
end
end
B=B/counter;
end
