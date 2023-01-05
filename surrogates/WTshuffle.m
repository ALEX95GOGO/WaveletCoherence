


%zero padding
exactlevels=log(length(x))/log(2);
numlevels=(floor(log(length(x))/log(2)));

count=1;
if abs(numlevels-exactlevels)>eps
    temp=zeros(2^(numlevels+1),1);
    count=2^(numlevels+1)-length(x)+1;
    temp(count:length(temp),1)=x;
    x=temp;
end

[cA,cD] = dwt(x,'mode',MODE)
