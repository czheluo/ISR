function []=ISR(famfile,genofile,outfile,n,p,ny,chr)
%default parameter
%if ~exist('outfile')
%    outfile = 'pop.traw.mat';
%end
%function ISR(arg1, arg2)
%end
%matlab -r "ISR('../data/pop.fam','../data/pop.traw','pop.traw.mat',87,28228,1,7);quit"
%example:ISR('../data/pop.fam','../data/pop.traw','pop.traw.mat',87,28228,1,7)
addpath([pwd,'/bin/']);
traw2mat(famfile,genofile,outfile,n,p,ny)
%traw2mat('../../pop.fam','../../poptraw.traw','popnew.mat',175,407045,3)
%pwd
%addpath('G:\D\D\MATLAB\isreg\rice2017\mouse\realanalysis\MLLM METHODS\ISR\bin')
%changing the names of output file with appropriate location
%resultname1='D:\\MATLAB\\isreg\\rice2017\\mouse\\realanalysis\\REAL_CFW\\RESULT\\adn_BMD(2chose)5.0.00001.txt';
%all_result_name='D:\\MATLAB\\isreg\\rice2017\\mouse\\realanalysis\\REAL_CFW\\RESULT\\all_BMD(2chose)5.0.00001.csv';
%opt_result_name='D:\\MATLAB\\isreg\\rice2017\\mouse\\realanalysis\\REAL_CFW\\RESULT\\aopt_BMD(2chose)5.0.00001.csv';
%genotype='D:/MATLAB/isreg/rice2017/thaliana/result/thaliana10snp6_genotype_file.xls';
%diary('ISR.log'); % notes the command window diary
%genotype file with completely numeric in .txt format 
%load mice_1161_97234.mat    
load(outfile);
%read the phenotype file with .csv or .xls Excel format
%kin=load('kinshipVan_lm.txt');%kin=kin(:,1:3);
%PCA=load('PCA5.txt');%PCA=PCA(:,1:3);
%impute=2;yname=2;
%y=xlsread('LD167.xlsx');%if the first of column and  
%y=xlsread('tha199_new_phe.xls');
%y=xlsread('mice_alldata.xls');
%yl=readtable('FileS2.txt');
%y=table2array(yl(:,2));
%y=xlsread('flc.xlsx');
%y=load('SDV.txt');


Y=y;
[~,p1]=size(Y);
if p1==1
    %clear,clc
    warning('off','all')
    diary('ISR.log'); % notes the command window diary
    impute=2;yname=2;
    alfa=0.00001;
    glm=5;%the number of running 
    ept0=2;%initially random chosed the number of multi-locus (probable result)
    gr=2;%export the result of significant association SNP genotype (with gr=1)
    nchr=chr;%the number of Chromosome
    sgt=0.05;%  bonferroni correction
    y=Y(:,1);
    %[xlm,textdata]=xlsread('rightdata/tha199_178384map.xlsx');
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
    ISR_REG
    %MLLM_COF
    %changing possible the last association results
    % writetable(opt_result,opt_result_name);writetable(all_result,all_result_name);
    diary off
    man=figure(1);
    set(man, 'Position', [200, 200, 1200, 450])
    manhattanhg;
    saveas(man,['manhattan.',num2str(p1),'.bmp']);
    clf;
    qq=figure(2);
    qqplot_conf;
    saveas(qq,['qq.',num2str(p1),'.bmp']);
    manqq=figure(3);
    subplot(1,2,1)
    manhattanhg
    subplot(1,2,2)
    qqplot_conf
    saveas(manqq,['manqq.',num2str(p1),'.bmp']);
elseif p1 > 1
    for i=1:p1
        impute=2;yname=2;
        diary(['isr',num2str(i),'.log']); % notes the command window diary
        alfa=0.00001;
        glm=5;%the number of running 
        ept0=2;%initially random chosed the number of multi-locus (probable result)
        gr=2;%export the result of significant association SNP genotype (with gr=1)
        nchr=chr;%the number of Chromosome
        sgt=0.05;%  bonferroni correction
        %[xlm,textdata]=xlsread('rightdata/tha199_178384map.xlsx');
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
        ISR_REG
        %MLLM_COF
        %changing possible the last association results
        % writetable(opt_result,opt_result_name);writetable(all_result,all_result_name);
        diary off
        man=figure(1);
        manhattanhg;
        saveas(man,['manhattan.',num2str(p1),'.bmp']);
        clf;
        qq=figure(2);
        qqplot_conf;
        saveas(qq,['qq.',num2str(p1),'.bmp']);
        manqq=figure(3);
        subplot(1,2,1)
        manhattanhg
        subplot(1,2,2)
        qqplot_conf
        saveas(manqq,['manqq.',num2str(p1),'.bmp']);
    end
end
end
