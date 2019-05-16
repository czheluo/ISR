
clear,clc,clf;
warning('off','all')

%genotype='D:/MATLAB/isreg/rice2017/thaliana/result/thaliana10snp6_genotype_file.xls';
%diary(resultname1) % notes the command window diary
%genotype file with completely numeric in .txt format 
%x=load('mice_1160_92734.txt');

%load t2d.mat
%label=labels;clear labels major minor
%x=full(X);clear X;
%read the phenotype file with .csv or .xls Excel format
%PCA=load('pca20.txt');PCA=PCA(:,1:5);npc=size(PCA,2);%PCA=PCA(:,1:3);
impute=2;yname=2;
%y=xlsread('tibial.xlsx');
%y=load('CFW_ALL_PHENO.txt');
%y=load('D1_2TOTDIST30.txt'); 
[y,labelna]=xlsread('D1_2TOTDIST30');t2=zeros(size(y,2),1);ylmn=y;
for ilmn=11:size(ylmn,2)
load cfw.mat
clf;
y=ylmn(:,ilmn);
%changing the names of output file with appropriate location
resultname1=(['tha/' 'adn_' char(labelna(1,ilmn)) '(5chose)3.0.00001.txt']);
all_result_name=(['tha/' 'all_' char(labelna(1,ilmn)) '(5chose)3.0.00001.csv']);
opt_result_name=(['tha/' 'opt_' char(labelna(1,ilmn)) '(5chose)3.0.00001.csv']);
%genotype='D:/MATLAB/isreg/rice2017/thaliana/result/thaliana10snp6_genotype_file.xls';
diary(resultname1) % notes the command window diary
%genotype file with completely numeric in .txt format 
%x=load('mice_1160_92734.txt');

%y=load('fastglucose.txt');
%y=y(:,3);
%choose the right model,1 or 2
mdl=1;

%ny=size(y,1);flm=zeros(ny,1);
%for i=1:ny
 %       if isnan(y(i,1))
  %          flm(i,1)=i;
   %     end
%end
%flm(flm==0)=[];x(flm,:)=[];y(isnan(y))=[];
glm=3;%the number of running 
ept0=5;%initially random chosed the number of multi-locus (probable result)
alfa=0.00001;%changing the singnificant level
gr=2;%export the result of significant association SNP genotype (with gr=1)
nchr=19;%the number of Chromosome
sgt=0.05;%  bonferroni correction
%[xlm,textdata]=xlsread('mousecfw_map.xlsx');
%label=textdata(:,1);chr=xlm(:,1);pos=xlm(:,3);
%label=textdata(1:1000,1);chr=xlm(1:1000,1);pos=xlm(1:1000,2);
%find the MISSING DATA add to the program
% [n,p]=size(x);
% for i=1:n
   %  for j=1:p
   %      if isnan(x(i,j))
  %          disp([i,j,x(i,j)])
  %      end
 %   end
% end
%stop
%SSR_rice_glmpc

%MLLM_Epistasis

MLLM
%MLLM_COF
%SSR_rice_P
%changing possible the last association results
%or write the txt file
writetable(opt_result,opt_result_name);writetable(all_result,all_result_name);

t2(i)=t;
clf;
man=figure(1);
Manhattanhg;
saveas(man,(['manhattanplot.' char(labelna(1,ilmn)) '(5chose)3_0.00001.bmp']));
clf;
qq=figure(2);
qqplot_conf;
saveas(qq,(['qqplot.' char(labelna(1,ilmn)) '(5chose)3_0.00001.bmp']));
diary off
%save the result in table in .txt
%[str,maxsize] = computer%To see the biggest memmory in your computer have
end

