Tout=LVTime(sout,cout,aout,T0,lastyear-firstyear-1);
yr=firstyear:1:lastyear;
figure
plot(yr,Tout,'.-','LineWidth',2);
hold
plot(yr(1:max(size(Data))), Data, 'sq','LineWidth',3)
grid
xlabel('Year')
ylabel('Units')
