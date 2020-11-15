function  GScluster(phefile,simufile, method, pdistArg)
%real phenotype cluster result
%GScluster('alltrait.xls','sence.xls', 'ward','mahalanobis')
%method linakge methods ward 
clear,clc;
figure(1)
[x,name]=xlsread(phefile);
cgo_all = clustergram(x','Colormap',parula,'RowPdist','chebychev','ColumnPdist', pdistArg,...
    'Linkage',method);%mahalanobis%chebychev
rn=name(2:20,1);cn=name(1,2:9);
set(cgo_all,'ColumnLabels',rn,'RowLabels',cn)
set(cgo_all,'Linkage','ward','Dendrogram',3)
set(cgo_all,'AnnotColor','k')
cgo_all.ColumnLabelsRotate=45;
figure(2)
%simulation cluster 
[x,name]=xlsread(simufile);
%cgo_all = clustergram(x','Colormap',parula,'RowPdist','chebychev','ColumnPdist','chebychev',...
  %  'Linkage','ward');%mahalanobis
cgo_all = clustergram(x','Colormap',parula,'ColumnPdist', pdistArg,...
    'Linkage',method);%mahalanobis%chebychev
rn=name(2:13,1);cn=name(1,2:9);
set(cgo_all,'ColumnLabels',rn,'RowLabels',cn)
set(cgo_all,'Linkage','ward','Dendrogram',3)
set(cgo_all,'AnnotColor','k')
cgo_all.ColumnLabelsRotate=45;
end










