clear,clc
warning('off','all')
addpath([pwd,'/src/']);
maxNumCompThreads(1);
load pop.mat
Y=y;X=x;
[n,p]=size(x);%rinf=zeros(nfolds,1);rref=zeros(nfolds,1);test=zeros(n,1);
ny=20;rinf=zeros(ny,1);rref=zeros(ny,1);%test=zeros(226,ny);
eff=zeros(p,ny);plmm=zeros(p,ny);tim=zeros(ny,1);beta=zeros(ny,1);
pve=zeros(p,ny);%ytest=zeros(226,ny);
for i=1:20
    t1=tic;
    N = size(X,1);
    sq=randperm(N);testing=sq(1:round(N*0.2));trainning=setdiff(sq,testing);
   %for i=1:nfolds
    %which=foldid==i;
    %if verbous, disp(['Fitting fold # ' num2str(i) ' of ' num2str(nfolds)]);end
   % x=reperma
    x=X(trainning,:);y=Y(trainning,1);
    [b,effect,plm,seblm,ft,lmbx,lm,r2lm,t] =ISR(x,y,2,2,0.01,3,25);
    x0=[ones(size(X(testing,:),1),1) X(testing,lm)];
    x1=[ones(size(X(trainning,:),1),1) X(trainning,lm)];
    %eff=zeros(size(x,2)+1,1);eff(1)=b(1); eff(2:end)=effect;
    predmat1=x0*b;ytest=predmat1;
    predmat2=x1*b;
    %rt=xo
    rinf(i)=corr(Y(testing,i),predmat1);
    rref(i)=corr(Y(trainning,i),predmat2);pve(:,i)=r2lm;
    %test(:,i)=testing;
    eff(:,i)=effect;plmm(:,i)=plm;tim(i)=toc(t1);beta(i)=b(1);
    save([num2str(i),'sw','.mat'],'rinf','rref','testing','eff','plmm','tim','beta','ytest','predmat2','t');
end
