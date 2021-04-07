h4=figure;
h4.Renderer='Painters';
subplot(2,1,1);
load kf.mat
manhattanplot (chr, pos, p)
subplot(2,1,2)
load yl.mat
manhattanplotinverted (chr, pos, p)
set(h4,'PaperSize',[10 10]); %set the paper size to what you want  
%print(h4,'filename','-dpdf')
print(h4,'filename','-dpdf','-fillpage')

