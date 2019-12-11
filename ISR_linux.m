function []=ISR_linux(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISR methods for GWAS @AUTHOR MENG LUO  %
% contact: czheluo@gmail.com             %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%seting the env path 
% alias matlab='/mnt/d/linux/MATLAB2016b/bin/matlab -nodesktop -nosplash -singleCompThread -logfile `date +%Y_%m_%d-%H_%M_%S`.log -r'
addpath([pwd,'/bin/']);
phefile = ft_getopt(varargin, 'phefile', 'phe.fam');
genofile = ft_getopt(varargin, 'genofile', 'pop.traw');
outfile = ft_getopt(varargin, 'outfile', 'pop.traw.mat');
sample = ft_getopt(varargin, 'sample',[]);
nSNP = ft_getopt(varargin, 'nSNP', []);
ntrait = ft_getopt(varargin, 'ntrait', []);
chr = ft_getopt(varargin, 'nchr', []);
opt_outresult = ft_getopt(varargin, 'opt_outresult', 'ISR.opt.outresult.txt');
all_outresult = ft_getopt(varargin, 'all_outresult', 'ISR.outresult.txt');
vcf = ft_getopt(varargin, 'vcf', []);
bed = ft_getopt(varargin, 'bed', []);
ncov = ft_getopt(varargin, 'ncov', []);
IM = ft_getopt(varargin, 'IM', 1);
sgv = ft_getopt(varargin, 'sgv', 0.05); % default bonferroni correction
% Usage:
%       matlab "ISR_linux('phefile','../data/pop.fam','genofile','../data/pop.traw','sample',87,'nSNP',28228,'ntrait',1,'ncov',5),exit;"
       
%default parameter
%if ~exist('outfile')
%    outfile = 'pop.traw.mat';
%end
%if ~exist('opt_outresult')
%    opt_outresult = 'opt_outresult.txt';
%end
%if ~exist('all_outresult')
%    all_outresult = 'all_outresult.txt';
%end
%maxNumCompThreads=1;
% write a shell script for ISR getopt....
if ~isempty(vcf)
    system(['plink ',' --vcf ',vcf,' --recode A-transpose --out pop']);
    traw2mat(phefile,genofile,outfile,sample,nSNP,ntrait,IM);
elseif ~isempty(bed)
    system(['plink ',' --bfile ',bed,' --recode A-transpose --out pop ']);
    traw2mat(phefile,genofile,outfile,sample,nSNP,ntrait,IM);
else
    traw2mat(phefile,genofile,outfile,sample,nSNP,ntrait,IM);
end
%traw2mat('../../pop.fam','../../poptraw.traw','popnew.mat',175,407045,3)
load(outfile);
Y=y;
[~,p1]=size(Y);
if p1==1
    %clear,clc
    warning('off','all')
    diary('ISR.log'); % notes the command window diary
    impute=2;yname=2;
    alfa=0.001;
    glm=3;%the number of running 
    ept0=5;%initially random chosed the number of multi-locus (probable result)
    gr=2;%export the result of significant association SNP genotype (with gr=1)
    nchr=chr;%the number of Chromosome
    sgt=sgv;%  bonferroni correction
    y=Y(:,p1);
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
    % ISR method
    if ~isempty(ncov)
        ISR_COV
    else
        ISR_REG
    end
    %changing possible the last association results
    writetable(opt_result,opt_outresult,'Delimiter','\t');writetable(all_result,all_outresult,'Delimiter','\t');
    diary off
    man=figure(1);
    set(man, 'Position', [200, 200, 1200, 450])
    manhattanhg;
    saveas(man,['manhattan.',num2str(p1),'.bmp']);
    clf;
    qq=figure(2);
    %set(qq, 'Position', [200, 200, 1200, 450])
    qqplot_conf;
    saveas(qq,['qq.',num2str(p1),'.bmp']);
    clf;
    %manqq=figure(3);
    %set(manqq, 'Position', [200, 200, 1200, 900])
    %subplot(2,1,1)
    %manhattanhg
    %subplot(2,1,2)
    %qqplot_conf
    %saveas(manqq,['manqq.',num2str(p1),'.bmp']);
    %clf;
elseif p1 > 1
    for i=1:p1
        impute=2;yname=2;
        diary(['isr',num2str(i),'.log']); % notes the command window diary
        alfa=0.00001;
        glm=5;%the number of running 
        ept0=2;%initially random chosed the number of multi-locus (probable result)
        gr=2;%export the result of significant association SNP genotype (with gr=1)
        nchr=chr;%the number of Chromosome
        sgt=sgv; 
        y=Y(:,p1);
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
        % ISR method
        if ~isempty(ncov)
            ISR_COV
        else
            ISR_REG
        end
        %changing possible the last association results
        writetable(opt_result,opt_outresult,'Delimiter','\t');writetable(all_result,all_outresult,'Delimiter','\t');
        diary off
        man=figure(1);
        set(man, 'Position', [200, 200, 1200, 450])
        manhattanhg;
        saveas(man,['manhattan.',num2str(p1),'.bmp']);
        clf;
        qq=figure(2);
        %set(qq, 'Position', [200, 200, 1200, 450])
        qqplot_conf;
        saveas(qq,['qq.',num2str(p1),'.bmp']);
        clf;
        %manqq=figure(3);
        %set(manqq, 'Position', [200, 200, 1200, 900])
        %subplot(2,1,1)
        %manhattanhg
        %subplot(2,1,2)
        %qqplot_conf
        %saveas(manqq,['manqq.',num2str(p1),'.bmp']);
        %clf;
    end
end
end
