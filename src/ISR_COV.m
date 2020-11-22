%ISR FOR GWAS with covariates
%@Author Meng Luo & Shiliang Gu  
%!copy fsdx.m fsd.m;
%deleted the missing phenotype (individuals)
%for l=n:-1:1
    %if isnan(y(l))
        %y(l)=[];
        %x(l,:)=[];
    %end
%end
ny=size(y,1);flm=zeros(ny,1);
display(['%number of total individuals= ',num2str(abs(ny))])
for i=1:ny
       if isnan(y(i,1))
           flm(i,1)=i;
      end
end
flm(flm==0)=[];x(flm,:)=[];y(isnan(y))=[];%PCA(flm,:)=[];
[n,p]=size(x);
%calculate the PCs as the Covariates
if ~isempty(ncov) && ncov > 1
    xy=zeros(p,1);
    for i = 1:p
        xy(i)=mean(x(:,i));
        x(:,i)=x(:,i) - xy(i);
    end
    [U,S,R]=svdk(x,ncov);
    pc=x*R;PCA=pc;npc=size(PCA,2);clear pc
end
display(['%number of analyzed individuals= ',num2str(abs(n))])%fprintf.....
display(['%number of analyzed SNPs= ',num2str(abs(p))])
display(['%number of covariates= ',num2str(abs(npc))])
plm=ones(p,1);effect=ones(p,1);r2lm=ones(p,1);seblm=ones(p,1);ft=ones(p,1);
%imputed genotype by mean
if impute==1
for j=1:p
    if find(isnan(xlz(:,j)))
        xnot=~isnan(xlz(:,j));xlzsum=sum(xnot);index=find(isnan(xlz(:,j)));
        xlz(index,j)=0;sumxlz=sum(xlz(:,j));mea=sumxlz/xlzsum;
        ave(1,j)=mea;med(1,j)=median(x(:,j));
    end
end 
 for i=1:n
     for j=1:p
         if isnan(x(i,j))
            fprintf(' Missing genotypes will be replaced with the mean genotype of SNP.\n')
            ylz=input('choose the impute methods = ')
            if ylz==1
                fprintf('imputed genotype by mean genotype of SNP.\n')
                x(i,j)=ave(j);
            elseif ylz==2
                     fprintf('imputed genotype by median genotype of SNP.\n')
                    x(i,j)=med(j);
            end
        end
    end
 end
end
%calculate the repeat phenotype observate value
%if yname(1,1)==yname(2,1) %the name with the cell name or char format
%[yn,yid,yidx]=unique(yname(:,1),'stable');
 %yval=accumarray(yidx,y(:,1),[],@mean); 
% y=accumarray(y(:,1), y(:,2), [], @mean);% if the name with numeric
 %y=yval;
%end
%[n,p]=size(x)
p1=p;
%alfa=.01;%change the alpha will change the power no significant
%mdl=input('Using Model II(without square term 2) or Model III(with square term 3) 2/3? ');
%mdl=1;
if mdl==1
    np=p1;
elseif mdl==2
    np=p1*(1+p1)/2;
elseif mdl==3
    np=p1*(3+p1)/2;
end
tr0=1:np;Tr0=tr0;tr1=1:p1;
XX=[ones(n,1),x];p0=p+1;
X=XX;
X0=zeros(n,np);
for i=1:p1
    X0(:,i)=x(:,i);
end
if mdl>1
    p2=p1;
    if mdl==3
        for i=1:p1
            p2=p2+1;
            X0(:,p2)=x(:,i).*x(:,i);
        end
    end
    for i=1:p1-1
        for j=i+1:p1
            p2=p2+1;
            X0(:,p2)=x(:,i).*x(:,j);
        end
    end
end
SSy=var(y)*(n-1);%alfa=.01;
% Stage I: primary regression
t1=tic;
ept=fix(1.5*ept0);
tr=randperm(p1);tr=tr(1:ept);%change the first chose
X=[ones(n,1),PCA,X0(:,tr)];p=size(X,2);%npc=size(PCA,2);
for lj=1:3
    dtr=setdiff(1:p1,tr); p3=length(dtr); cp=randperm(p3);
    for li=1:p3
        p=p+1;
        X(:,p)=X0(:,dtr(cp(li)));tr(p-npc-1)=dtr(cp(li));
        if rank(X)<p% change the p number 
            X(:,p)=[];tr(p-npc-1)=[];p=p-1;
        end
        a=X'*X;c=inv(a);k=X'*y;b=a\k;
        q=y'*y-b'*k;mse=q/(n-p);
        up=b.*b./diag(c);up(1)=[]; r2=(SSy-q)/SSy;
        f=up/mse;flq=f(npc+1:end);pr=fcdf(min(flq),1,n-p,'upper');
        qi=find(flq-min(flq)==0);
        %calculate the p value 
        SEb=sqrt(mse*diag(c));SEb(1)=[];seblm(tr,1)=SEb(npc+1:end);
        flz=f(npc+1:end);uplz=up(npc+1:end);
        plz=fcdf(flz,1,n-p,'upper');plm(tr)=plz;ft(tr)=flz;
        r2lz=uplz./SSy;effect(tr)=b(npc+2:end);r2lm(tr)=r2lz;%calculate the PVE
        if pr>=falf1(alfa/p,p+fr21(r2)) && p>npc+5%chose more SNPs with simulation
            tr(qi)=[];X(:,qi+npc+1)=[];p=p-1;
        end
    end
end
a=X'*X;c=inv(a);k=X'*y;b=a\k;q=y'*y-b'*k;mse=q/(n-p);
up=b.*b./diag(c);up(1)=[];r2=(SSy-q)/SSy;
f=up/mse;
%calculate the p value 
flz=f(npc+1:end);uplz=up(npc+1:end);pr=fcdf(min(flz),1,n-p,'upper');
plz=fcdf(flz,1,n-p,'upper');plm(tr)=plz;ft(tr)=flz;
r2lz=uplz./SSy;effect(tr)=b(npc+2:end);r2lm(tr)=r2lz;%calculate the PVE
disp(['Initial regression, ','p=',num2str(p-npc-1),', R^2=',num2str(r2)])
% Stage II: Re-examining all effect terms repeatedly!
dtr=setdiff(Tr0,tr);dtr1=[];sg=-1;fj=0;pt=0;
alf1=.8;cr4=1;alf=.1;rp=0;rp1=2;rp2=1;v=-1;Of=-100;OF=-1000;
nc=0;ne=0;cr=.25;ecr=.15;btr=[];etr=[];nj=0;FR=rand(1,16)-1000;TrR=zeros(16,50);
while rp<=glm %can change
    p3=length(dtr); cp=randperm(p3);
    ii=0;v=v+1;if v>13+.5*rp,v=1;end
    pct3=.15+pt+.001*mod(rp,150); pct4=.85-pt-.001*mod(rp,150);
    while ii<p3
        ii=ii+1;
        flq=f(npc+1:end);qi=find(flq-min(flq)==0);%find the nsigificant tr
        if p>2.7*ept+mod(v,27)+.03*mod(rp,350)*(randn+.1)
            if pr<alf1,dtr1=union(dtr1,tr(qi));end
            tr(qi)=[];X(:,qi+npc+1)=[];p=p-1;
        elseif p<=5+.01*p1
            nc=nc+1;
        else
            alf=cr*alfa/(1+.25*p+.25*fr21(r2));
            alf=.001*(alf<.001)+(alf>=.001)*alf*(alf<=.6)+.6*(alf>.6);
            if pr>=alf
                if pr<alf1,dtr1=union(dtr1,tr(qi));end
                tr(qi)=[];X(:,qi+npc+1)=[];p=p-1;
            else
                nc=nc+1;
            end
        end
        p=p+1;X(:,p)=X0(:,dtr(cp(ii)));tr(p-npc-1)=dtr(cp(ii));
        if rank(X)<p
            X(:,p)=[];tr(p-npc-1)=[];p=p-1;
        end
        a=X'*X;c=inv(a);k=X'*y;
        b=a\k;q=y'*y-b'*k;mse=q/(n-p);
        up=b.*b./diag(c);up(1)=[];
        f=up/mse; SEblz=sqrt(mse*diag(c));SEblz(1)=[];seblm(tr,1)=SEblz(npc+1:end);
        %calculate the p value
        flz=f(npc+1:end);uplz=up(npc+1:end);pr=fcdf(min(flz),1,n-p,'upper');
        plz=fcdf(flz,1,n-p,'upper');plm(tr)=plz;ft(tr)=flz;
        pl=1+length(tr)+.5*sum(tr>p1);
        r2=(SSy-q)/SSy;of=fsd(pl/n)*r2/(p-1)/(1-r2)*(n-p);
        r2lz=uplz./SSy;effect(tr)=b(npc+2:end);r2lm(tr)=r2lz;%calculate the PVE
        if of>Of
            Bx=X;Btr=tr;Of=of;btr=union(btr,tr);
            fof=str2double(num2str(of,7));
            if sum(FR==fof)==0
                FR(16)=fof; TrR(16,1:end)=0; TrR(16,1:length(Btr))=Btr;
                [FR,ind]=sort(FR,'descend');
                TrR=TrR(ind,:);
            end
        elseif of>Of-ecr
            ne=ne+1;etr=union(etr,tr);
        end  
        alfa=0.01;
        if length(f)==npc%pc proplems
           fprintf('Not Have variable to choose.\n');
           break
        else
        while p>2.6*ept+mod(v,25) || pr>alfa+.018*mod(v,63)
            if pr<falf1(alf1/p,p+fr21(r2)) || p<=1.25*ept+.5*mod(v,22)+.01*mod(rp,190)*randn
                break
            end
             flq=f(npc+1:end);qj=find(flq-min(flq)==0);
            rd=rand(1);
            if rd<=.2+pt && rp>=5 && p-npc-1>=10
                [F,ind]=sort(flq,'descend');
                re=randperm(round((p-npc-1)*pct3));q=round((p-npc-1)*pct4-.05)+re(1);
                qi=ind(q);
            else
                qi=qj;
            end
            if pr<alf1,dtr1=union(dtr1,tr(qi));end
            X(:,qi+npc+1)=[];p=p-1;tr(qi)=[];
            a=X'*X;c=inv(a);k=X'*y;b=a\k;
            q=y'*y-b'*k;mse=q/(n-p);
            up=b.*b./diag(c);up(1)=[];
            f=up/mse;SEblz=sqrt(mse*diag(c));SEblz(1)=[];seblm(tr,1)=SEblz(npc+1:end);
           %calculate the p value 
            flz=f(npc+1:end);uplz=up(npc+1:end);pr=fcdf(min(flz),1,n-p,'upper');
            plz=fcdf(flz,1,n-p,'upper');plm(tr)=plz;ft(tr)=flz;%calculate the p value
            pl=1+length(tr)+.5*sum(tr>p1);
            r2=(SSy-q)/SSy;of=fsd(pl/n)*r2/(p-1)/(1-r2)*(n-p);
            r2lz=uplz./SSy;effect(tr)=b(npc+2:end);r2lm(tr)=r2lz;%calculate the PVE
            if of>Of
                Bx=X;Btr=tr;Of=of;btr=union(btr,tr);
                fof=str2double(num2str(of,7));
                if sum(FR==fof)==0
                    FR(16)=fof; TrR(16,1:end)=0; TrR(16,1:length(Btr))=Btr;
                    [FR,ind]=sort(FR,'descend');
                    TrR=TrR(ind,:);
                end
            elseif of>Of-ecr
                ne=ne+1;etr=union(etr,tr);
            end
        end
        end
    end
    p5=p;
    pause(.0001)
    pct1=.18+pt+.001*mod(rp+20,220); pct2=.82-pt-.001*mod(rp+20,220);
    if mod(rp,135)>75, rd1=random('unif',.58,1.18,1); else rd1=.68; end
    alfa=0.001;
    while p>rd1*ept+.5*sg*mod(v,20)+.025*mod(rp,180)*randn || pr>=alfa
        if p<=3+.005*p1
            break
        elseif pr<falf1(alf1/p,p+fr21(r2)) && p<=.25*ept+.25*mod(v,17)
            break
        end
        flq=f(npc+1:end);qj=find(flq-min(flq)==0);
        rd=rand(1);
        if rd<=.25+pt && rp>5
            [F,ind]=sort(flq,'descend');
            re=randperm(round((p-npc-1)*pct1));
            q=round((p-npc-1)*pct2-.05)+re(1);
            qi=ind(q);
        else
            qi=qj;
        end
        X(:,qi+1)=[];dtr1=union(dtr1,tr(qi));tr(qi)=[];p=p-1;
        a=X'*X;c=inv(a);k=X'*y;b=a\k;
        q=y'*y-b'*k;mse=q/(n-p);
        up=b.*b./diag(c);up(1)=[];
        f=up/mse;SEblz=sqrt(mse*diag(c));SEblz(1)=[];seblm(tr,1)=SEblz(npc+1:end);
        %calculate the p value
        flz=f(npc+1:end);uplz=up(npc+1:end);pr=fcdf(min(flz),1,n-p,'upper');
        plz=fcdf(flz,1,n-p,'upper');plm(tr)=plz;ft(tr)=flz;
        pl=1+length(tr)+.5*sum(tr>p1);
        r2=(SSy-q)/SSy;of=fsd(pl/n)*r2/(p-1)/(1-r2)*(n-p);
        r2lz=uplz./SSy;effect(tr)=b(npc+2:end);r2lm(tr)=r2lz;%calculate the PVE
        if of>Of
            Bx=X;Btr=tr;Of=of;btr=union(btr,tr);
            fof=str2double(num2str(of,7));
            if sum(FR==fof)==0
                FR(16)=fof; TrR(16,1:end)=0; TrR(16,1:length(Btr))=Btr;
                [FR,ind]=sort(FR,'descend');
                TrR=TrR(ind,:);
            end
        elseif of>Of-ecr
            ne=ne+1;etr=union(etr,tr);
        end
    end
    p6=p;
    pause(.001)
    for lj=1:2
        rd=rand(1);
        if rd<.6
            Dtr=union(btr,1:p1);
        elseif rd>.85
            if length(etr)>.125*np, etr=btr;end
            Dtr=union(etr,1:p1);
        else
            Dtr=1:p1;
        end
        dtr=setdiff(Dtr,tr); p3=length(dtr); cp=randperm(p3);
        for li=1:p3
            p=p+1;
            X(:,p)=X0(:,dtr(cp(li)));tr(p-npc-1)=dtr(cp(li));
            if rank(X)<p
                dtr1=union(dtr1,tr(p-npc-1));
                X(:,p)=[];tr(p-npc-1)=[];p=p-1;
            end
            a=X'*X;c=inv(a);k=X'*y;b=a\k;
            q=y'*y-b'*k;mse=q/(n-p);
            up=b.*b./diag(c);up(1)=[];
            f=up/mse;SEblz=sqrt(mse*diag(c));SEblz(1)=[];seblm(tr,1)=SEblz(npc+1:end);
            %calculate the p value
            flz=f(npc+1:end);uplz=up(npc+1:end);pr=fcdf(min(flz),1,n-p,'upper');
            plz=fcdf(flz,1,n-p,'upper');plm(tr)=plz;ft(tr)=flz;
            pl=1+length(tr)+.5*sum(tr>p1);
            r2=(SSy-q)/SSy;of=fsd(pl/n)*r2/(p-1)/(1-r2)*(n-p);
            r2lz=uplz./SSy;effect(tr)=b(npc+2:end);r2lm(tr)=r2lz;%calculate the PVE
            if of>Of
                Bx=X;Btr=tr;Of=of;btr=union(btr,tr);
                fof=str2double(num2str(of,7));
                if sum(FR==fof)==0
                    FR(16)=fof; TrR(16,1:end)=0; TrR(16,1:length(Btr))=Btr;
                    [FR,ind]=sort(FR,'descend');
                    TrR=TrR(ind,:);
                end
            elseif of>Of-ecr/100
                ne=ne+1;etr=union(etr,tr);
            end
            flq=f(npc+1:end);qi=find(flq-min(flq)==0);%change  the f value 
            if pr>=5*falf1(alfa/p,p+fr21(r2)) && p>5
                dtr1=union(dtr1,tr(qi));
                tr(qi)=[];X(:,qi+npc+1)=[];p=p-1;
                a=X'*X;k=X'*y;b=a\k;c=inv(a);
                q=y'*y-b'*k;
                %add by Meng
                mse=q/(n-p);up=b.*b./diag(c);up(1)=[];f=up/mse;
                pl=1+length(tr)+.5*sum(tr>p1);
                r2=(SSy-q)/SSy;of=fsd(pl/n)*r2/(p-1)/(1-r2)*(n-p);
                %calculate the p value 
                SEblz=sqrt(mse*diag(c));SEblz(1)=[];seblm(tr,1)=SEblz(npc+1:end);
               flz=f(npc+1:end);uplz=up(npc+1:end);%pr=1-fcdf(min(flz),1,n-p);%change the pr crit only the SNPs
               plz=fcdf(flz,1,n-p,'upper');plm(tr)=plz;blz=b(npc+2:end);ft(tr)=flz;
               r2lz=uplz./SSy;effect(tr)=blz;r2lm(tr)=r2lz;%calculate the PVE
                if of>Of
                    Bx=X;Btr=tr;Of=of;btr=union(btr,tr);
                    fof=str2double(num2str(of,7));
                    if sum(FR==fof)==0
                        FR(16)=fof; TrR(16,1:end)=0; TrR(16,1:length(Btr))=Btr;
                        [FR,ind]=sort(FR,'descend');
                        TrR=TrR(ind,:);
                    end
                elseif of>Of-ecr/100
                    ne=ne+1;etr=union(etr,tr);
                end
            end
        end
    end
    %p7=p;disp([p5,p6,p7])
    pause(.001), nj=nj+1;
    if v<=0
        btr=[];etr=[];
    elseif length(btr)<=1
        btr=Btr;btr(btr>p1)=[];
    end
    if mod(v,2)==1
        cr1=exp(-6.3*(nc-10-.01*p3)/(.1*p3+.1*p1));
        cr1=.1*(cr1<.1)+(cr1>=.1)*cr1*(cr1<=10)+10*(cr1>10);
        cr=.25*cr*(1+cr1)+.25*(cr+cr1);
        cr2=.65*exp(-5.2*(ne-8-.0065*p3)/(.065*p3+.1*p1))+.3*exp(-.7*(length(etr)-2.5*ept-.1*np./(1+12.5*exp(-.025*nj)))/(.01*np+2.5*nj));
        cr2=.1*(cr2<.1)+(cr2>=.1)*cr2*(cr2<=8)+8*(cr2>8);
        ecr=.24*ecr*(1+cr2)+.24*(ecr+cr2);
    else
        cr1=exp(-6.1*(nc-18-.02*p3)/(.1*p3+.1*p1));
        cr1=.1*(cr1<.1)+(cr1>=.1)*cr1*(cr1<=10)+10*(cr1>10);
        cr=.25*cr*(1+cr1)+.25*(cr+cr1);
        cr2=.65*exp(-5*(ne-15-.013*p3)/(.065*p3+.1*p1))+.3*exp(-.7*(length(etr)-2.5*ept-.1*np./(1+12.5*exp(-.025*nj)))/(.01*np+2.5*nj));
        cr2=.1*(cr2<.1)+(cr2>=.1)*cr2*(cr2<=8)+8*(cr2>8);
        ecr=.24*ecr*(1+cr2)+.24*(ecr+cr2);
        nc=0;ne=0;
    end
    if mod(v,15)==5+fj && rp>=2 % using current result restart searching all terms
        dtr=setdiff(Tr0,tr); p3=length(dtr)+.15;
        alf1=falf2(p3); dtr1=[];
    elseif mod(v,15)==10+fj && rp>=2 % using previous best result restart searching all terms
        r3=randperm(3);tr3=TrR(r3(1),:); tr4=TrR(r3(2),:); rd=rand;
        if rd<.33, tr=union(tr3,tr4); elseif rd>.67, tr=intersect(tr3,tr4); else tr=Btr; end
        tr(tr==0)=[];X=[ones(n,1),PCA,X0(:,tr)]; p=size(X,2);dtr=setdiff(Tr0,tr);
        dtr1=[]; p3=length(dtr); alf1=falf2(p3)+.15;ecr=.65*ecr;
    elseif mod(v,15)==0 && rp>=7
        if mod(rp1,7)==1 % intensive search for terms appeared in good results
            rd=rand;
            if rd<.33, tr=union(tr,Btr); elseif rd>.67 tr=intersect(tr,Btr);else tr=Btr;end
            X=[ones(n,1),PCA,X0(:,tr)];p=size(X,2);
            Dtr=union(btr,etr);dtr=setdiff(Dtr,tr); p3=length(dtr);
            alf1=falf2(p3)+.25;dtr1=[];
            rp1=rp1+1;if length(etr)>100+.035*np;etr=btr;nj=5;end
            sg=-1;fj=0;pt=0;
            disp(['searching 1: intensive search based on Btr, p=',num2str(p-npc-1)])
        elseif mod(rp1,7)==2 %restart search all terms based on part terms appeared in btr
            p44=fix(1.25*ept+.01*mod(rp,330)*rand);
            while p44>.85*length(btr), p44=p44-1; end
            X=[ones(n,1),PCA,X0(:,btr)]; tr=btr; p=size(X,2);X=X+randn(n,p)/10000;
            X(:,1)=ones(n,1);
            while p>p44
                a=X'*X;c=inv(a);k=X'*y;b=a\k; q=y'*y-b'*k;mse=q/(n-p);
                up=b.*b./diag(c);up(1)=[]; f=up/mse;
                flq=f(npc+1:end);qi=find(flq-min(flq)==0);
                X(:,qi+npc+1)=[];tr(qi)=[];p=p-1;
            end
            X=[ones(n,1),PCA,X0(:,tr)];p=size(X,2);
            dtr=setdiff(Tr0,tr); p3=length(dtr);alf1=falf2(p3)+.15;fj=0;
            dtr1=[];rp1=rp1+1;pt=0;
            rd=rand;if rd<.25, sg=-.5; else sg=-1; end
            disp(['searching 2: restart search based on partial(stepwise) btr, p=',num2str(p44+1-npc)]);
        elseif mod(rp1,7)==3 % intensive search based on terms appeared in best results
            tr3=unique(TrR(1:7,:)); tr4=tr3'; tr4(tr4==0)=[]; p4=length(tr4);
            p44=p4-1; re=randperm(p4); tr=tr4(re(1:p44));
            X=[ones(n,1),PCA,X0(:,tr)];p=size(X,2);
            Dtr=union(btr,etr);Dtr=union(Dtr,1:p1);
            dtr=setdiff(Dtr,tr); p3=length(dtr);alf1=falf2(p3)+.25;
            dtr1=[];rp1=rp1+1;fj=0;pt=0;
            rd=rand;if rd<.25, sg=-.25; else sg=-.5;end
            if length(etr)>200+.1*np;etr=btr;nj=5;end
            disp(['searching 3: intensive search based on top 4 Btr, p=',num2str(p-npc-1)]);
        elseif mod(rp1,7)==4
            if length(Btr)<.85*ept,ce=1.15;elseif length(Btr)>1.15*ept,ce=.85;else ce=1.05;end
            r3=randperm(3); tr2=TrR(r3(1),:); tr3=TrR(r3(2),:);
            l1=length(intersect(tr,intersect(tr2,tr3)));
            l2=length(intersect(tr,union(tr2,tr3)));
            l3=length(union(tr,intersect(tr2,tr3)));
            l4=length(union(tr,union(tr2,tr3)));
            d=[l1,l2,l3,l4];e=ce*ept+.03*randn(1)*mod(rp,110);
            se=find(abs(d-e)-min(abs(d-e))==0);
            if se==1
                tr=intersect(tr,intersect(tr2,tr3));
            elseif se==2
                tr=intersect(tr,union(tr2,tr3));
            elseif se==3
                tr=union(tr,intersect(tr2,tr3));
            else
                tr=union(tr,union(tr2,tr3));
            end
            tr(tr==0)=[];
            X=[ones(n,1),PCA,X0(:,tr)];p=size(X,2);dtr=setdiff(Tr0,tr);
            p3=length(dtr);alf1=falf2(p3)+.05;sg=rand-1;
            dtr1=[];rp1=rp1+1; fj=0;pt=.05;
            disp(['searching 4: restart search based on top 5 results, p=',num2str(p-npc-1)]);
        elseif mod(rp1,7)==5 % restart search based on partial etr terms
            tr3=unique(TrR(1:9,:)); tr4=tr3'; tr4(tr4==0)=[]; p4=length(tr4);
            if p4>10, p44=p4-2; else p44=p4-1;end
            re=randperm(p4); tr=tr4(re(1:p44));
            X=[ones(n,1),PCA,X0(:,tr)];p=size(X,2);
            Dtr=union(btr,etr);Dtr=union(Dtr,1:p1);
            dtr=setdiff(Dtr,tr); p3=length(dtr);alf1=falf2(p3)+.75;
            dtr1=[];rp1=rp1+1; sg=-1;fj=15;pt=.15;
            if length(etr)>150+.067*np;etr=btr;nj=5;end
            disp(['searching 5: intensive search based on top 6 Btrs, p=',num2str(p-npc-1)]);
        elseif mod(rp1,7)==6
            p44=fix(1.5*ept+.015*mod(rp,150)*rand);
            while p44>.85*length(btr),p44=p44-1;end
            p4=length(btr);re=randperm(p4);
            tr=btr(re(1:p44)); 
            X=[ones(n,1),PCA,X0(:,tr)]; p=size(X,2);fj=15;pt=.15;
            dtr=setdiff(Tr0,tr); p3=length(dtr);alf1=falf2(p3)+.25;
            dtr1=[];rp1=rp1+1;sg=-1;
            if mod(rp2,2)==1
                Of=(.68+.05*randn)*Of;ecr=.15*ecr;rp2=rp2+1;
            else
                Of=(.5+.05*randn)*Of;ecr=.1*ecr;btr=[];etr=[];nj=0;rp2=rp2+1;
            end
            disp(['searching 6: restart search based on partial etr, p=',num2str(p-npc-1)]);
        else  %restart search all terms based on partial best terms
            tr3=unique(TrR); tr4=tr3'; tr4(tr4==0)=[]; tr=tr4;
            p44=fix(.8*ept)+randi(5,1); 
            while p44>.8*length(tr),p44=p44-1;end
            X=[ones(n,1),PCA,X0(:,tr)]; p=size(X,2);
            while p>p44
                a=X'*X;c=inv(a);k=X'*y;b=a\k; q=y'*y-b'*k;mse=q/(n-p);
                up=b.*b./diag(c);up(1)=[]; f=up/mse;
                SEblz=sqrt(mse*diag(c));SEblz(1)=[];seblm(tr,1)=SEblz(npc+1:end);
                %calculate the p value
                flz=f(npc+1:end);uplz=up(npc+1:end);
                plz=fcdf(flz,1,n-p,'upper');plm(tr)=plz;ft(tr)=flz;
                pl=1+length(tr)+.5*sum(tr>p1);
                r2=(SSy-q)/SSy;of=fsd(pl/n)*r2/(p-1)/(1-r2)*(n-p);
                r2lz=uplz./SSy;effect(tr)=b(npc+2:end);r2lm(tr)=r2lz;%calculate the PVE
                if of>Of
                    Bx=X;Btr=tr;Of=of;btr=union(btr,tr);
                    fof=str2double(num2str(of,7));
                    if sum(FR==fof)==0
                        FR(16)=fof; TrR(16,1:end)=0; TrR(16,1:length(Btr))=Btr;
                        [FR,ind]=sort(FR,'descend');
                        TrR=TrR(ind,:);
                    end
                end
                flq=f(npc+1:end);qi=find(flq-min(flq)==0);
                X(:,qi+npc+1)=[];tr(qi)=[];p=p-1;
            end
            disp(['searching 7: restart search based on tr4(-3), p=',num2str(p44+1-npc)]);
            X=[ones(n,1),PCA,X0(:,tr)]; p=size(X,2);fj=15;pt=.15;
            dtr=setdiff(Tr0,tr); p3=length(dtr);alf1=falf2(p3)+.3*rand;
            dtr1=[];rp1=rp1+1;sg=-1;
            Of=(.88+.01*randn)*Of;ecr=.2*ecr;
        end
    else % intensive search for the weak effect terms
        dtr=setdiff(dtr1,tr);
        p3=length(dtr); alf1=falf2(p3)+.15;
        dtr1=[];
    end
    while rank(X)<p
        dtr1=union(dtr1,tr(p-npc-1));
        X(:,p)=[];tr(p-npc-1)=[];p=p-1;
    end
    a=X'*X;c=inv(a);k=X'*y;b=a\k;
    q=y'*y-b'*k;mse=q/(n-p);
    up=b.*b./diag(c);up(1)=[];
    f=up/mse;
    SEblz=sqrt(mse*diag(c));SEblz(1)=[];seblm(tr,1)=SEblz(npc+1:end);
    %calculate the p value 
    flz=f(npc+1:end);uplz=up(npc+1:end);pr=fcdf(min(flz),1,n-p,'upper');
    plz=fcdf(flz,1,n-p,'upper');plm(tr)=plz;ft(tr)=flz;
    r2=(SSy-q)/SSy;r2lz=uplz./SSy;effect(tr)=b(npc+2:end);r2lm(tr)=r2lz;%calculate the PVE
    pl=1+length(tr)+.5*sum(tr>p1);
    of=fsd(pl/n)*r2/(p-1)/(1-r2)*(n-p);
    if of>Of
        Bx=X;Btr=tr;Of=of;btr=union(btr,tr);
        fof=str2double(num2str(of,7));
        if sum(FR==fof)==0
            FR(16)=fof; TrR(16,1:end)=0; TrR(16,1:length(Btr))=Btr;
            [FR,ind]=sort(FR,'descend');
            TrR=TrR(ind,:);
        end
    elseif of>Of-ecr
        ne=ne+1;etr=union(etr,tr);
    end
    if Of>OF
        OF=Of;rp=0;BX=Bx;BTR=Btr;
    else
        rp=rp+1;
    end
    ept=fix(.8*ept0+.065*mod(rp,165));
    qb=y'*y-(Bx\y)'*Bx'*y;r2=(SSy-qb)/SSy;p4=size(Bx,2);p20=length(btr);p30=length(etr);
    disp(['v=',num2str(v),', p=',num2str(p4-npc-1),', R^2=',num2str(r2),', of=',num2str(Of),',  btrsize=',num2str(p20),', etrsize=',num2str(p30),', ecr=',num2str(ecr)])
end
X=BX;tr=BTR;p=size(X,2);
a=X'*X;k=X'*y;c=inv(a);b=a\k;
q=y'*y-b'*k;mse=q/(n-p);
up=b.*b./diag(c);up(1)=[];
f=up/mse;SEblz=sqrt(mse*diag(c));SEblz(1)=[];seblm(tr,1)=SEblz(npc+1:end);
%calculate the p value 
flz=f(npc+1:end);uplz=up(npc+1:end);pr=fcdf(flz,1,n-p,'upper');
plz=fcdf(flz,1,n-p,'upper');plm(tr)=plz;ft(tr)=flz;%calculate the p value
SEb=sqrt(mse*diag(c));SEb(1)=[];
r2lz=uplz./SSy;effect(tr)=b(npc+2:end);r2lm(tr)=r2lz;%calculate the PVE
%Stage IV display last results.
t=toc(t1);
if mdl==1
        fo=0;
        %p1=17000
end
for i=1:p1
    if i<10
        TR(i,:)=char(['X',num2str(i),'     ']);
    elseif i>=10 && i<100
        TR(i,:)=char(['X',num2str(i),'    ']);
    elseif i>=100 && i<1000
        TR(i,:)=char(['X',num2str(i),'   ']);
    elseif i>=1000 && i<10000
        TR(i,:)=char(['X',num2str(i),'  ']);
    elseif i>=10000 && i<100000
        TR(i,:)=char(['X',num2str(i),' ']);
    elseif i>=100000 && i<1000000
        TR(i,:)=char(['X',num2str(i)]);
    else
        TR(i,:)=char(['X',num2str(i)]);
    end
end
Tr=TR(tr,:);
disp(['    p=',num2str(p-npc-1)])
disp('Last Results:')
disp('Xi       bi       SEb       Up        F       p')
disp(['X0      ',num2str(b(1))])
for i=1:p-npc-1
    disp([Tr(i,:),'  ',num2str(b(i+1)),'  ',num2str(SEb(i)),'  ',num2str(up(i)),'  ',num2str(f(i)),'  ',num2str(pr(i))])
end
disp(['Error  ',num2str(n-p),' ',num2str(q),' ',num2str(mse)])
disp(['Total  ',num2str(n-1),' ' num2str(SSy)])
r2=(SSy-q)/SSy;disp(['p=',num2str(p-npc-1),', R^2=',num2str(r2),', OF=',num2str(OF)])
dfe=n-p;dfT=n-1;
r2lz=uplz./SSy;effect(tr)=b(npc+2:end);r2lm(tr)=r2lz;%calculate the PVE
[n1,m1]=size(Tr);
for i=1:n1
 new_tr(i,:)=strrep(Tr(i,:),'X',' ');
end
lm=str2num(new_tr);markers=marker(lm);chrgl=chr(lm);poslm=pos(lm);
Marker=markers;chromosome=chrgl;position=poslm;p_f=pr;beta=b(npc+2:end);sebeta=SEb(npc+1:end);r2=r2lz;
opt_result=table(Marker,chromosome,position,flz,p_f,beta,sebeta,r2);
SNP=marker;Chromosome=chr;Position=pos;Ft=ft;P_F=plm;Beta=effect;SEbeta=seblm;R2=r2lm;
all_result=table(SNP,Chromosome,Position,Ft,P_F,Beta,SEbeta,R2);
save('ISR_result.mat','all_result','opt_result','marker','chr','pos','plm','-v7.3');

lmbx=x(:,lm);
if gr==1
xlswrite(genotype,lmbx,1);
%xlswrite(genotype,Bx,2);
end
