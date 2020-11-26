filename = '.\M6K4_2.txt';
[e1,e2,e3,e4]=textread(filename,'%f%f%f%f','delimiter',',');

Y= [mean(e1);mean(e2);mean(e3);mean(e4)];
b=diag(Y);
c=bar(b,0.5,'stack');
color=[1,0,1;0.62745,0.12549,0.94118;1,0.64706,0;0.80392,0.78824,0.78824];
for i=1:4
set(c(i),'FaceColor',color(i,:));
end
Y_1 = Y;
for i = 1:length(Y)
    text(i,Y_1(i),num2str(Y_1(i)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'FontName','Times New Roman');
end
ylabel('Localization Error (km)');
xlabel('Standard Deviation of TDOA Measurement Noise (1000 ns)');
set(gca,'xtick',1:4);
set(gca,'XTickLabel',{'Emitter 1','Emitter 2','Emitter 3','Emitter 4'},'FontSize',12,'FontName','Times New Roman'); 