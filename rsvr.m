function errs=rsvr(n,epi)
u=rand(n,1)*(1+epi);
ry=[];
for i=1:n
    if u(i)<epi
        ry(i)=u(i);
    else
        ry(i)=epi-log(1+epi-u(i));
    end
end
    errs=ry.*(2*randi(2,1,n)-3);
    errs=errs';
end

