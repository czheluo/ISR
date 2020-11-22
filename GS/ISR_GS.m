clear,clc
warning('off','all') 
load rice-HDRA0.1.mat x y
load rice_maf3.mat 
mafn=find(maf3<0.01);
 Y=y;X=x(:,mafn);
[n,p]=size(x);%rinf=zeros(nfolds,1);rref=zeros(nfolds,1);test=zeros(n,1);
ny=1;rinf=zeros(ny,1);rref=zeros(ny,1);test=zeros(226,ny);
eff=zeros(p,ny);plmm=zeros(p,ny);tim=zeros(ny,1);beta=zeros(ny,1);
pve=zeros(p,ny);%ytest=zeros(226,ny);
for ny=1:10
    t=tic;
     N = size(X,1);
    sq=randperm(N);testing=sq(1:round(N*0.2));trainning=setdiff(sq,testing);
   %for i=1:nfolds
    %which=foldid==i;
    %if verbous, disp(['Fitting fold # ' num2str(i) ' of ' num2str(nfolds)]);end
   % x=reperma
    x=X(trainning,:);y=Y(trainning,ny);
    [b,effect,plm,seblm,ft,lmbx,lm,r2lm] =ISR(x,y,2,2,0.01,3,25);
    %Y,x,y,impute,yname,alfa,glm,ept0
    %cvfit =glmnet(x(~which,:), y(~which),family, options);
    %predmat(which,:) = glmnetPredict(cvfit, type,x(which,:),options.lambda);
    x0=[ones(size(X(testing,:),1),1) X(testing,lm)];
    x1=[ones(size(X(trainning,:),1),1) X(trainning,lm)];
    %eff=zeros(size(x,2)+1,1);eff(1)=b(1); eff(2:end)=effect;
    predmat1=x0*b;ytest=predmat1;
    predmat2=x1*b;
    %rt=xo
    rinf(ny)=corr(Y(testing,ny),predmat1);
    rref(ny)=corr(Y(trainning,ny),predmat2);pve(:,ny)=r2lm;
    test(:,ny)=testing;eff(:,ny)=effect;plmm(:,ny)=plm;tim(ny)=toc(t);beta(ny)=b(1);
    save('rice_3_25_root_(1_10)_120922.mat','rinf','rref','test','eff','plmm','tim','beta','ytest');
end
%end

