function D = makesurftemp(Elev)

% create look up table for elevation dependent temperature pofile for
% ice sheet temperatures

load('tempvselev.mat') % 'the file with the temperatures and elevations'

goodInd = ~isnan(Elev);
D = polyfit(Elev(goodInd),Temp1(goodInd),2);

figure(1)
plot(Elev,Temp1,'.','Markersize',10)

x = linspace(min(Elev),max(Elev),50);
f1 = polyval(D,x);

figure(1)
hold on
plot(x,f1,'linewidth',3)
xlabel('Elevation above sea level (m)')
ylabel('Temperature (^oC)')

print('Surfacetemp','-dpng')

end



