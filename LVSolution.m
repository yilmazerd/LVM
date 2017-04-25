[p k] = size(Data);

a_lowest = rand(p,1);
s_lowest = ones(p,p) - rand(p,p)*2;
c_lowest =rand(p,p)/10;
sum_lowest = 1e9;


for i = 1:iterations,
a = rand(p,1);
s = ones(p,p) - rand(p,p)*2;
c =rand(p,p)/10;
[sout,cout,aout]=LVoptimize(Data,true,a,s,c);


T0=Data(:,1);
Tout=LVTime(sout,cout,aout,T0,lastyear-firstyear-1);
yr=firstyear:1:lastyear;



disp('Total sum of squares')
sumsq = sum(sum((Data(1:end,2:end) - Tout(1:p,2:k)).^2))/((p-1)*k);

if (sumsq<sum_lowest)
    sum_lowest = sumsq;
    a_lowest = aout;
    s_lowest = sout;
    c_lowest = cout;
end

end

aout = a_lowest;
sout = s_lowest;
cout = c_lowest;