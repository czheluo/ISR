%@MENG LUO
%CONTACT chzeluo@gmail.com
clear,clc;
load na23_chr4_2.mat
%tit='Chromosome 4';
subplot(3,1,1)
bin=30000;%windon size 
p1=max(po);p2=min(po);
p=round((p1-p2)/bin);pl=zeros(length(po),1);
%for i=1:1
%a=0;
pol=po;
pm=zeros(p,1);
for j=1:p
    [a,b]=find(po<(po(1)+bin));
    pm(j)=length(a);
    col=zeros(p,3);
        %pl(a>pm(j))=j;%
    cl=unifrnd(0,1,1,3);
    col(j,:)=cl;
    stem(po(a),.05*ones(size(po(a))),'Color',cl,'Marker','none');
        %a=a+0;
    po(a)=[];  
    hold on
end
pm=num2str(pm);
lgd = legend(pm,'Location','north','Orientation','horizontal');
title(lgd,'Plotted SNPs');
%h2=stem(po,.05*ones(size(po)),'Color','k','Marker','none');
set(gca,'Xtick',[],'Ytick',[]);ylim([0 0.05]);%ylabel('plotted SNPs')
set(gca,'FontName','Times New Roman','FontWeight','bold','FontSize',14);

subplot(3,1,2)
pml=-log10(pml);cb=30;
scatter(pol,pml,cb,pml,'filled','LineWidth',2)
%colormap jet% the catgory of colormap{parula,jet,hsv,hot,cool,spring,summer,autumn,winter,gray,bone,copper,pink,lines,colorcube,prism,flag,white}
cm=colorbar('location','east');
cm.Label.String = 'r^2';
%new color type
colormap(cbrewer2('RdBu'));
box on;
set(gca,'FontName','Times New Roman', 'FontWeight','bold','FontSize',14);
ylabel('-log_{10}(\itP)');%title(tit);
set(gca,'XTick',[]);

subplot(3,1,3);
for i=1:16
    a1=[al(i),al(i)];
    ps=px(i,:)/1000000;
    %line(ps,a1,'LineWidth',7);
    %annotation('arrow',)
    %dim=[ps,a1];
    %annotation(hf, 'arrow', dim);
    text(mean(ps)-0.005,a1(1)+0.3,lab(i));
    %davinci( 'arrow', 'X', ps, 'Y',a1,'ArrowType', 'double','LineWidth',2);%,'Shaft.Type','rectangle'
    %daspect( [1 1 1] ) 
    davinci( 'arrow', 'X', ps, ...
                  'Y',a1, ...
                  'ArrowType','single',...%double
                  'Shaft.Type','rectangle', ...
                  'Shaft.Width',0.3,... #'Head.Length',0.03,...
                  'Head.Width',0.3, ...
                  'Head.Sweep',0, ...
                  'Color','b', ...
                  'EdgeColor','b', ...
                  'LineWidth',0.5 );
    hold on
    %setting the right or left arrow (if 0 or 1)
    %text(ps,a1,'\rightarrow','FontSize',14,'FontWeight','bold')
end
set(gca,'FontName','Times New Roman', 'FontWeight','bold','FontSize',14);
%set(gca,'xlim',[pos(n1),pos(n2)]);
set(gca,'YTick',[]);
xlabel('Position on Chr(Mb)');set(gca,'Box','on','XGrid','on');   





