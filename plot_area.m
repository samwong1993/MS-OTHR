clear all
index = [0 100 200 300 400 500 600 700 800 900 1000];
for i = 1:11
    filename = ".\M6K4_"+string(index(i))+".txt";
    [e1,e2,e3,e4]=textread(filename,'%f%f%f%f','delimiter',',');
    Y(i,:) = [mean(e1),mean(e2),mean(e3),mean(e4)];
end
h = area(index,max(Y'),"facecolor","r","edgecolor",[1,1,1])
hold on
set(h,'edgealpha',0,'facealpha',0.3)
area(index,min(Y'),"facecolor",[1,1,1],"edgecolor",[1,1,1])
% plot(mean(Y'),'b')
h1 = plot(index,mean(Y'),'ok-','linewidth',1.1,'markerfacecolor',[36, 169, 225]/255);
legend([h,h1],"Error Range","Mean Error")
xlabel('Standard Deviation of TDOA Measurement Noise (ns)')
ylabel('RMSE of the Geolocation(m)')

