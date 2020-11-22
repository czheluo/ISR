%clear,clc
%@author Meng Luo qqplot with confidence interval 
%p=csvread('book4.csv');p=p(1:2710,1);
%plm=csvread('P_D2hact0to30.csv');
%plm=csvread('p.csv');
o = -log10(sort(plm,'ascend'));%c=6;o=o(1:100000);
alphaLevel=0.05;halflevel=0.5;oneMinusAlphalevel=1-alphaLevel;%95%confidience
x=1:1:length(o);
e = -log10( (x-0.5)./(length(o)))';
cl=load('colorchrhg.txt');
nclm=randperm(length(cl));
scatter(e,o,30,cl(nclm(1),:),'o','filled')%6,8,12,15,17,20,
hold on
M=length(o);
mSeq=10.^(log10(0.5):0.1:(log10(M - 0.5)+0.1));
n=length(mSeq);n1=length(mSeq);alpha=ones(n,1);half=ones(n,1);oneMinusAlpha=ones(n,1);
if alphaLevel==alphaLevel
    for i=1:n
        alpha(i)=betainv(alphaLevel,mSeq(i),M-mSeq(i));
    end
end
if halflevel==halflevel
    for i=1:n
        half(i)=betainv(halflevel,mSeq(i),M-mSeq(i));
    end
end
if oneMinusAlphalevel==1-alphaLevel
    for i=1:n
        oneMinusAlpha(i)=betainv(1-alphaLevel,mSeq(i),M-mSeq(i));
    end
end
%alpha=[alpha half oneMinusAlpha];
betaDown=half-alpha;betaUp=oneMinusAlpha-half; 
theoreticalPval=(mSeq/M)';
lowerBar=-log10(theoreticalPval-betaDown);%lowerBar=real(lowerBar);
upperBar=-log10(theoreticalPval+betaUp);
yBar=-log10(theoreticalPval);
%pn=theoreticalPval-betaDown;
%ny=length(yBar);
%for i=1:ny
    %if pn(i)<0
       % lowerBar(i)=[];upperBar(i)=[];yBar(i)=[];
   % end
%end       
%lowerBar=real(lowerBar);
g=[0 0.5 0];
line(yBar,lowerBar,'Color',g,'LineStyle','-.','LineWidth',2.5)
line(yBar,upperBar,'Color',g,'LineStyle','-.','LineWidth',2.5)
line(e,e,'Color','k','LineStyle','-','LineWidth',2.5)
set(gca,'FontName','Times New Roman','FontSize',16);
xlabel('Expected -log_{10}(\itP)','FontName','Times New Roman','FontSize',16);
ylabel('Observed -log_{10}(\itP)','FontName','Times New Roman','FontSize',16);
legend('ISR','Location','northwest')
%logQuantiles=-log10(e);
%N=length(p);
%n1 = length(o);
%c95=zeros(n1,1);c05=zeros(n1,1);
%for i=1:n1          
        %    xi = ceil((10^-logQuantiles(i)) * N);
            %c95(i,1) =betainv(0.95, xi, N - xi + 1);
          %  c05(i,1)=betainv(0.05, xi, N - xi + 1);
        %    if xi == 0
        %      xi = 1;
         %      c95(i,1) =betainv(0.95, xi, N - xi + 1);
        %       c05(i,1)=betainv(0.05, xi, N - xi + 1);
       %     end
%end
%hold on
%x=0:0.1:5;
%plot(e,c95)
