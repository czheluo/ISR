%MLLM FOR GWAS
%@Author Shiliang Gu & Meng Luo
%!copy fsdx.m fsd.m;
%[n,p]=size(x)
%[n,p]=size(x);
%display('An Efficient Repeatedly Screening (Stepwise Screening) Multi-locus Linear Model Approach For GWAS')
%display(['%number of analyzed individuals= ',num2str(abs(n))])%fprintf.....
%display(['%number of analyzed SNPs= ',num2str(abs(p))])
%NOTE: Now we may take only the significant SNPs according to their marginal epistatic effects and run a simple exhaustive search between them
ps=find(plm<(0.05/(size(x,2))));
%xps=x(:,ps);
%plm=ones(p,1);effect=ones(p,1);r2lm=ones(p,1);seblm=ones(p,1);ft=ones(p,1);
%p1=p;%alfa=.01;%change the alpha will change the power no significant
%mdl=input('Using Model II(without square term 2) or Model III(with square term 3) 2/3? ');
%mdl=1;
np=length(ps)*(size(x,2)-1);
if np==0
    fprintf('no singnificant varians detected.\n')
    return
end
display(['%number of analyzed effects = ',num2str(abs(np))])
%scan for dominance or additive epistatic effects
trnum=zeros(np,2);Eplm=ones(np,1);Eeffect=ones(np,1);
Er2lm=ones(np,1);Eseblm=ones(np,1);Eft=ones(np,1);
fo=0;
if mdl==2
    %for i=1:p
       % fo=fo+1;
       % trnum(fo,1)=i;trnum(fo,2)=i;%single locus model
    %end
    for i=1:length(ps)
        for j=1:size(x,2)
            fo=fo+1;%change the position to start search
                trnum(fo,1)=ps(i);trnum(fo,2)=j;
        end
    end
    for it=1:length(trnum)
        if trnum(it,1)==trnum(it,2)
            itr(it)=it;
        else
            itr(it)=0;
        end        
    end
    itr(itr==0)=[];
    trnum(itr,:)=[];
end
tr0=1:np;Tr0=tr0;tr1=1:p1;
%XX=[ones(n,1),x];
p0=p+1;
%X=XX;
%X0=zeros(n,np);
%scan for dominance or additive epistatic effects
%for i=1:p1
  %  X0(:,i)=x(:,i);
%end
%if mdl==2
    %scan for dominance-by-dominance or additive-by-additive epistatic effects
  %  p2=0;
    %if mdl==3
        %for i=1:p1
            %p2=p2+1;
            %X0(:,p2)=x(:,i).*x(:,i);
        %end
    %end
   % for i=1:p1-1
    %    for j=i+1:p1
     %       p2=p2+1;
     %       X0(:,p2)=x(:,i).*x(:,j);
    %    end
  %  end
%else
    % scan for dominance-by-additive epistatic effects
    %p2=p1/2;p3=0;
    %for i=1:p1-1
     %   for j=p2+1:p1
      %      p3=p3+1;
       %     X0(:,p3)=x(:,i).*x(:,j);
       % end
   % end
%end
%calculated the consuming time
t1=tic;
SSy=var(y)*(n-1);%alfa=.01;
% Stage I: primary regression
ept=fix(1.5*ept0);
tr=randperm(length(trnum));tr=tr(1:ept);%change the first chose
for i=1:length(tr)
    %X0=zeros(n,length(tr));
        X0(:,i)=x(:,trnum(tr(i),1)).*x(:,trnum(tr(i),2));
end
X=[ones(n,1),X0];p=size(X,2);p1=length(trnum);
for lj=1:3
    dtr=setdiff(1:p1,tr); p3=length(dtr); cp=randperm(p3);
    for li=1:p3
        p=p+1;
        %for i=1:length(dtr(cp(li)))
            %X0=zeros(n,length(tr)):
            X0=x(:,trnum(dtr(cp(li)),1)).*x(:,trnum(dtr(cp(li)),2));
        %end
        %X0=x(:,trnum(dtr(cp(li)),1)).*x(:,trnum(dtr(cp(li)),2));
        X(:,p)=X0;tr(p-1)=dtr(cp(li));
        if rank(X)<p
            X(:,p)=[];tr(p-1)=[];p=p-1;
        end
        a=X'*X;c=inv(a);k=X'*y;b=a\k;
        q=y'*y-b'*k;mse=q/(n-p);
        up=b.*b./diag(c);up(1)=[]; Er2=(SSy-q)/SSy;
        f=up/mse;pr=fcdf(min(f),1,n-p,'upper');
        SEb=sqrt(mse*diag(c));SEb(1)=[];Eseblm(tr,1)=SEb;
        qi=find(f-min(f)==0);
        %calculate the p value 
       %nlz=size(PCA,2);
        plz=fcdf(f,1,n-p,'upper');Eplm(tr)=plz;Eft(tr)=f;
        r2lz=up./SSy;Eeffect(tr)=b(2:end);Er2lm(tr)=r2lz;%calculate the PVE
        if pr>=falf1(alfa/p,p+fr21(Er2)) && p>5%choose more SNPs here to change
            tr(qi)=[];X(:,qi+1)=[];p=p-1;
        end
    end
end
a=X'*X;c=inv(a);k=X'*y;b=a\k;q=y'*y-b'*k;mse=q/(n-p);
up=b.*b./diag(c);up(1)=[];Er2=(SSy-q)/SSy;
f=up/mse;pr=fcdf(min(f),1,n-p,'upper');
SEblz=sqrt(mse*diag(c));SEblz(1)=[];Eseblm(tr,1)=SEblz;
%calculate the p value 
plz=fcdf(f,1,n-p,'upper');Eplm(tr)=plz;Eft(tr)=f;
r2lz=up./SSy;Eeffect(tr)=b(2:end);Er2lm(tr)=r2lz;%calculate the PVE
disp(['Initial regression, ','p=',num2str(p),', R^2=',num2str(Er2)])
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
        qi=find(f-min(f)==0);
        if p>2.7*ept+mod(v,27)+.03*mod(rp,350)*(randn+.1)
            if pr<alf1,dtr1=union(dtr1,tr(qi));end
            tr(qi)=[];X(:,qi+1)=[];p=p-1;
        elseif p<=5+.01*p1
            nc=nc+1;
        else
            alf=cr*alfa/(1+.25*p+.25*fr21(Er2));
            alf=.001*(alf<.001)+(alf>=.001)*alf*(alf<=.6)+.6*(alf>.6);
            if pr>=alf
                if pr<alf1,dtr1=union(dtr1,tr(qi));end
                tr(qi)=[];X(:,qi+1)=[];p=p-1;
            else
                nc=nc+1;
            end
        end
        X0=x(:,trnum(dtr(cp(ii)),1)).*x(:,trnum(dtr(cp(ii)),2));
        p=p+1;X(:,p)=X0;tr(p-1)=dtr(cp(ii));
        if rank(X)<p
            X(:,p)=[];tr(p-1)=[];p=p-1;
        end
        a=X'*X;c=inv(a);k=X'*y;
        b=a\k;q=y'*y-b'*k;mse=q/(n-p);
        up=b.*b./diag(c);up(1)=[];
        f=up/mse; pr=fcdf(min(f),1,n-p,'upper');
        SEblz=sqrt(mse*diag(c));SEblz(1)=[];Eseblm(tr,1)=SEblz;
        plz=fcdf(f,1,n-p,'upper');Eplm(tr)=plz;Eft(tr)=f;%calculate the p value
        pl=1+length(tr)+.5*sum(tr>p1);
        Er2=(SSy-q)/SSy;of=fsd(pl/n)*Er2/(p-1)/(1-Er2)*(n-p);
        r2lz=up./SSy;Eeffect(tr)=b(2:end);Er2lm(tr)=r2lz;%calculate the PVE
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
        while p>2.6*ept+mod(v,25) || pr>alfa+.018*mod(v,63)
            if pr<falf1(alf1/p,p+fr21(Er2)) || p<=1.25*ept+.5*mod(v,22)+.01*mod(rp,190)*randn
                break
            end
            qj=find(f-min(f)==0);
            rd=rand(1);
            if rd<=.2+pt && rp>=5
                [F,ind]=sort(f,'descend');
                re=randperm(round((p-1)*pct3));q=round((p-1)*pct4-.05)+re(1);
                qi=ind(q);
            else
                qi=qj;
            end
            if pr<alf1,dtr1=union(dtr1,tr(qi));end
            X(:,qi+1)=[];p=p-1;tr(qi)=[];
            a=X'*X;c=inv(a);k=X'*y;b=a\k;
            q=y'*y-b'*k;mse=q/(n-p);
            up=b.*b./diag(c);up(1)=[];
            f=up/mse;pr=fcdf(min(f),1,n-p,'upper');
            SEblz=sqrt(mse*diag(c));SEblz(1)=[];Eseblm(tr,1)=SEblz;
            plz=fcdf(f,1,n-p,'upper');Eplm(tr)=plz;Eft(tr)=f;%calculate the p value
            pl=1+length(tr)+.5*sum(tr>p1);
            Er2=(SSy-q)/SSy;of=fsd(pl/n)*Er2/(p-1)/(1-Er2)*(n-p);
            r2lz=up./SSy;Eeffect(tr)=b(2:end);Er2lm(tr)=r2lz;%calculate the PVE
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
    p5=p;
    pause(.0001)
    pct1=.18+pt+.001*mod(rp+20,220); pct2=.82-pt-.001*mod(rp+20,220);
    if mod(rp,135)>75, rd1=random('unif',.58,1.18,1); else rd1=.68; end
    alfa=0.001;
    while p>rd1*ept+.5*sg*mod(v,20)+.025*mod(rp,180)*randn || pr>=alfa
        if p<=3+.005*p1
            break
        elseif pr<falf1(alf1/p,p+fr21(Er2)) && p<=.25*ept+.25*mod(v,17)
            break
        end
        qj=find(f-min(f)==0);
        rd=rand(1);
        if rd<=.25+pt && rp>5
            [F,ind]=sort(f,'descend');
            re=randperm(round((p-1)*pct1));
            q=round((p-1)*pct2-.05)+re(1);
            qi=ind(q);
        else
            qi=qj;
        end
        X(:,qi+1)=[];dtr1=union(dtr1,tr(qi));tr(qi)=[];p=p-1;
        a=X'*X;c=inv(a);k=X'*y;b=a\k;
        q=y'*y-b'*k;mse=q/(n-p);
        up=b.*b./diag(c);up(1)=[];
        f=up/mse;pr=fcdf(min(f),1,n-p,'upper');
        SEblz=sqrt(mse*diag(c));SEblz(1)=[];Eseblm(tr,1)=SEblz;
        plz=fcdf(f,1,n-p,'upper');Eplm(tr)=plz;Eft(tr)=f;%calculate the p value
        pl=1+length(tr)+.5*sum(tr>p1);
        Er2=(SSy-q)/SSy;of=fsd(pl/n)*Er2/(p-1)/(1-Er2)*(n-p);
        r2lz=up./SSy;Eeffect(tr)=b(2:end);Er2lm(tr)=r2lz;%calculate the PVE
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
            X0=x(:,trnum(dtr(cp(li)),1)).*x(:,trnum(dtr(cp(li)),2));
            X(:,p)=X0;tr(p-1)=dtr(cp(li));
            if rank(X)<p
                dtr1=union(dtr1,tr(p-1));
                X(:,p)=[];tr(p-1)=[];p=p-1;
            end
            a=X'*X;c=inv(a);k=X'*y;b=a\k;
            q=y'*y-b'*k;mse=q/(n-p);
            up=b.*b./diag(c);up(1)=[];
            f=up/mse;pr=fcdf(min(f),1,n-p,'upper');
            SEblz=sqrt(mse*diag(c));SEblz(1)=[];Eseblm(tr,1)=SEblz;
            plz=fcdf(f,1,n-p,'upper');Eplm(tr)=plz;Eft(tr)=f;%calculate the p value
            pl=1+length(tr)+.5*sum(tr>p1);
            Er2=(SSy-q)/SSy;of=fsd(pl/n)*Er2/(p-1)/(1-Er2)*(n-p);
            r2lz=up./SSy;Eeffect(tr)=b(2:end);Er2lm(tr)=r2lz;%calculate the PVE
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
            qi=find(f-min(f)==0);
            if pr>=5*falf1(alfa/p,p+fr21(Er2)) && p>5
                dtr1=union(dtr1,tr(qi));
                tr(qi)=[];X(:,qi+1)=[];p=p-1;
                a=X'*X;k=X'*y;b=a\k;
                q=y'*y-b'*k;
                pl=1+length(tr)+.5*sum(tr>p1);
                Er2=(SSy-q)/SSy;of=fsd(pl/n)*Er2/(p-1)/(1-Er2)*(n-p);
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
        if rd<.33, tr=union(tr3,tr4);tr(tr==0)=[]; elseif rd>.67, tr=intersect(tr3,tr4);tr(tr==0)=[]; else tr=Btr;tr(tr==0)=[]; end
        for ijk=1:length(tr)
            %X0=size(n,length(tr));
                X0(:,ijk)=x(:,trnum(tr(ijk),1)).*x(:,trnum(tr(ijk),2));
        end
        tr(tr==0)=[];X=[ones(n,1),X0]; p=size(X,2);dtr=setdiff(Tr0,tr);
        dtr1=[]; p3=length(dtr); alf1=falf2(p3)+.15;ecr=.65*ecr;
    elseif mod(v,15)==0 && rp>=7
        if mod(rp1,7)==1 % intensive search for terms appeared in good results
            rd=rand;
            if rd<.33, tr=union(tr,Btr); elseif rd>.67 tr=intersect(tr,Btr);else tr=Btr;end
            for ijk=1:length(tr)
                X0(:,ijk)=x(:,trnum(tr(ijk),1)).*x(:,trnum(tr(ijk),2));
            end
            %X0=x(:,trnum(tr,1)).*x(:,trnum(tr,2));
            X=[ones(n,1),X0];p=size(X,2);
            Dtr=union(btr,etr);dtr=setdiff(Dtr,tr); p3=length(dtr);
            alf1=falf2(p3)+.25;dtr1=[];
            rp1=rp1+1;if length(etr)>100+.035*np;etr=btr;nj=5;end
            sg=-1;fj=0;pt=0;
            disp(['searching 1: intensive search based on Btr, p=',num2str(p)])
        elseif mod(rp1,7)==2 %restart search all terms based on part terms appeared in btr
            p44=fix(1.25*ept+.01*mod(rp,330)*rand);
            while p44>.85*length(btr), p44=p44-1; end
            for ijk=1:length(btr)
                X0(:,ijk)=x(:,trnum(btr(ijk),1)).*x(:,trnum(btr(ijk),2));
            end
            %X0=x(:,trnum(btr,1)).*x(:,trnum(btr,2));
            X=[ones(n,1),X0]; tr=btr; p=size(X,2);X=X+randn(n,p)/10000;
            X(:,1)=ones(n,1);
            while p>p44
                a=X'*X;c=inv(a);k=X'*y;b=a\k; q=y'*y-b'*k;mse=q/(n-p);
                up=b.*b./diag(c);up(1)=[]; f=up/mse;
                qi=find(f-min(f)==0);
                X(:,qi+1)=[];tr(qi)=[];p=p-1;
            end
            for ijk=1:length(tr)
                X0(:,ijk)=x(:,trnum(tr(ijk),1)).*x(:,trnum(tr(ijk),2));
            end
            %X0=x(:,trnum(tr,1)).*x(:,trnum(tr,2));                        
            X=[ones(n,1),X0];p=size(X,2);
            dtr=setdiff(Tr0,tr); p3=length(dtr);alf1=falf2(p3)+.15;fj=0;
            dtr1=[];rp1=rp1+1;pt=0;
            rd=rand;if rd<.25, sg=-.5; else sg=-1; end
            disp(['searching 2: restart search based on partial(stepwise) btr, p=',num2str(p44+1)]);
        elseif mod(rp1,7)==3 % intensive search based on terms appeared in best results
            tr3=unique(TrR(1:7,:)); tr4=tr3'; tr4(tr4==0)=[]; p4=length(tr4);
            p44=p4-1; re=randperm(p4); tr=tr4(re(1:p44));
            for ijk=1:length(tr)
                X0(:,ijk)=x(:,trnum(tr(ijk),1)).*x(:,trnum(tr(ijk),2));
            end
            %X0=x(:,trnum(tr,1)).*x(:,trnum(tr,2));
            X=[ones(n,1),X0];p=size(X,2);
            Dtr=union(btr,etr);Dtr=union(Dtr,1:p1);
            dtr=setdiff(Dtr,tr); p3=length(dtr);alf1=falf2(p3)+.25;
            dtr1=[];rp1=rp1+1;fj=0;pt=0;
            rd=rand;if rd<.25, sg=-.25; else sg=-.5;end
            if length(etr)>200+.1*np;etr=btr;nj=5;end
            disp(['searching 3: intensive search based on top 4 Btr, p=',num2str(p)]);
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
            for ijk=1:length(tr)
                X0(:,ijk)=x(:,trnum(tr(ijk),1)).*x(:,trnum(tr(ijk),2));
            end
            %X0=x(:,trnum(tr,1)).*x(:,trnum(tr,2));
            X=[ones(n,1),X0];p=size(X,2);dtr=setdiff(Tr0,tr);
            p3=length(dtr);alf1=falf2(p3)+.05;sg=rand-1;
            dtr1=[];rp1=rp1+1; fj=0;pt=.05;
            disp(['searching 4: restart search based on top 5 results, p=',num2str(p)]);
        elseif mod(rp1,7)==5 % restart search based on partial etr terms
            tr3=unique(TrR(1:9,:)); tr4=tr3'; tr4(tr4==0)=[]; p4=length(tr4);
            if p4>10, p44=p4-2; else p44=p4-1;end
            re=randperm(p4); tr=tr4(re(1:p44));
            for ijk=1:length(tr)
                X0(:,ijk)=x(:,trnum(tr(ijk),1)).*x(:,trnum(tr(ijk),2));
            end
            %X0=x(:,trnum(tr,1)).*x(:,trnum(tr,2));
            X=[ones(n,1),X0];p=size(X,2);
            Dtr=union(btr,etr);Dtr=union(Dtr,1:p1);
            dtr=setdiff(Dtr,tr); p3=length(dtr);alf1=falf2(p3)+.75;
            dtr1=[];rp1=rp1+1; sg=-1;fj=15;pt=.15;
            if length(etr)>150+.067*np;etr=btr;nj=5;end
            disp(['searching 5: intensive search based on top 6 Btrs, p=',num2str(p)]);
        elseif mod(rp1,7)==6
            p44=fix(1.5*ept+.015*mod(rp,150)*rand);
            while p44>.85*length(btr),p44=p44-1;end
            p4=length(btr);re=randperm(p4);
            tr=btr(re(1:p44)); 
            for ijk=1:length(tr)
                X0(:,ijk)=x(:,trnum(tr(ijk),1)).*x(:,trnum(tr(ijk),2));
            end
            %X0=x(:,trnum(tr,1)).*x(:,trnum(tr,2));
            X=[ones(n,1),X0]; p=size(X,2);fj=15;pt=.15;
            dtr=setdiff(Tr0,tr); p3=length(dtr);alf1=falf2(p3)+.25;
            dtr1=[];rp1=rp1+1;sg=-1;
            if mod(rp2,2)==1
                Of=(.68+.05*randn)*Of;ecr=.15*ecr;rp2=rp2+1;
            else
                Of=(.5+.05*randn)*Of;ecr=.1*ecr;btr=[];etr=[];nj=0;rp2=rp2+1;
            end
            disp(['searching 6: restart search based on partial etr, p=',num2str(p)]);
        else  %restart search all terms based on partial best terms
            tr3=unique(TrR); tr4=tr3'; tr4(tr4==0)=[]; tr=tr4;
            p44=fix(.8*ept)+randi(5,1); 
            while p44>.8*length(tr),p44=p44-1;end
            for ijk=1:length(tr)
                X0(:,ijk)=x(:,trnum(tr(ijk),1)).*x(:,trnum(tr(ijk),2));
            end
            %X0=x(:,trnum(tr,1)).*x(:,trnum(tr,2));     
            X=[ones(n,1),X0]; p=size(X,2);
            while p>p44
                a=X'*X;c=inv(a);k=X'*y;b=a\k; q=y'*y-b'*k;mse=q/(n-p);
                up=b.*b./diag(c);up(1)=[]; f=up/mse;
                SEblz=sqrt(mse*diag(c));SEblz(1)=[];Eseblm(tr,1)=SEblz;
                plz=fcdf(f,1,n-p,'upper');Eplm(tr)=plz;Eft(tr)=f;%calculate the p value
                pl=1+length(tr)+.5*sum(tr>p1);
                Er2=(SSy-q)/SSy;of=fsd(pl/n)*Er2/(p-1)/(1-Er2)*(n-p);
                r2lz=up./SSy;Eeffect(tr)=b(2:end);Er2lm(tr)=r2lz;%calculate the PVE
                if of>Of
                    Bx=X;Btr=tr;Of=of;btr=union(btr,tr);
                    fof=str2double(num2str(of,7));
                    if sum(FR==fof)==0
                        FR(16)=fof; TrR(16,1:end)=0; TrR(16,1:length(Btr))=Btr;
                        [FR,ind]=sort(FR,'descend');
                        TrR=TrR(ind,:);
                    end
                end
                qi=find(f-min(f)==0);
                X(:,qi+1)=[];tr(qi)=[];p=p-1;
            end
            disp(['searching 7: restart search based on tr4(-3), p=',num2str(p44+1)]);
            for ijk=1:length(tr)
                X0(:,ijk)=x(:,trnum(tr(ijk),1)).*x(:,trnum(tr(ijk),2));
            end
            %X0=x(:,trnum(tr,1)).*x(:,trnum(tr,2));     
            X=[ones(n,1),X0]; p=size(X,2);fj=15;pt=.15;
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
        dtr1=union(dtr1,tr(p-1));
        X(:,p)=[];tr(p-1)=[];p=p-1;
    end
    a=X'*X;c=inv(a);k=X'*y;b=a\k;
    q=y'*y-b'*k;mse=q/(n-p);
    up=b.*b./diag(c);up(1)=[];
    f=up/mse;pr=fcdf(min(f),1,n-p,'upper');
    SEblz=sqrt(mse*diag(c));SEblz(1)=[];Eseblm(tr,1)=SEblz;
    plz=fcdf(f,1,n-p,'upper');Eplm(tr)=plz;Eft(tr)=f;%calculate the p value
    Er2=(SSy-q)/SSy;r2lz=up./SSy;Eeffect(tr)=b(2:end);Er2lm(tr)=r2lz;%calculate the PVE
    pl=1+length(tr)+.5*sum(tr>p1);
    of=fsd(pl/n)*Er2/(p-1)/(1-Er2)*(n-p);
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
    qb=y'*y-(Bx\y)'*Bx'*y;Er2=(SSy-qb)/SSy;p4=size(Bx,2);p20=length(btr);p30=length(etr);
    disp(['v=',num2str(v),', p=',num2str(p4),', R^2=',num2str(Er2),', of=',num2str(Of),',  btrsize=',num2str(p20),', etrsize=',num2str(p30),', ecr=',num2str(ecr)])
end
X=BX;tr=BTR;p=size(X,2);
a=X'*X;k=X'*y;c=inv(a);b=a\k;
q=y'*y-b'*k;mse=q/(n-p);
up=b.*b./diag(c);up(1)=[];
f=up/mse;pr=fcdf(f,1,n-p,'upper');
SEblz=sqrt(mse*diag(c));SEblz(1)=[];Eseblm(tr,1)=SEblz;
plz=fcdf(f,1,n-p,'upper');Eplm(tr)=plz;Eft(tr)=f;%calculate the p value
SEb=sqrt(mse*diag(c));SEb(1)=[];
r2lz=up./SSy;Eeffect(tr)=b(2:end);Er2lm(tr)=r2lz;%calculate the PVE
t=toc(t1);

%Stage IV display last results.

%npr=size(pr,1);
%ppr=sort(pr);ppr(ppr==0)=[];
%for i=1:npr
%if pr(i,1)==0
  %  pr(i,:)=min(ppr);
%end
%end
%Tr=TR(tr,:);

disp(['    p=',num2str(p)])
disp('Last Results:')
disp(['b0        ',num2str(b(1))])
disp('SNP1      SNP2       bi       SEb       Up        F       p')
for i=1:p-1
    disp([['SNP' num2str(trnum(tr(i),1))],'     ',['SNP' num2str(trnum(tr(i),2))],'    ',num2str(b(i+1)),'  ',num2str(SEb(i)),'  ',num2str(up(i)),'  ',num2str(f(i)),'  ',num2str(pr(i))])
end
disp(['    ']);
disp(['Error  ',num2str(n-p),' ',num2str(q),' ',num2str(mse)]);
disp(['Total  ',num2str(n-1),' ' num2str(SSy)]);
disp(['    ']);
Er2=(SSy-q)/SSy;disp(['p=',num2str(p),', R^2=',num2str(Er2),', OF=',num2str(OF)]);
dfe=n-p;dfT=n-1;
r2lz=up./SSy;Eeffect(tr)=b(2:end);Er2lm(tr)=r2lz;%calculate the PVE
nplm=size(Eplm,1);
%for i=1:nplm
%if plm(i,1)==0
  %  plm(i,2)=fpdf(f(i,2),1,n-p);
    %plm(i,1)=min(pr);
%end
%end  
%nplm=size(plm,1);
%for i=1:nplm
  %  if plm(i,1)==plm(tr,1)
       % plm(i,1)=plm(tr,1);
   % elseif plm(i,1)~=plm(tr,1) %&& plm(i,1)<0.05/nplm
    %    plm(i,1)=0;
   % end
%end
%[n1,m1]=size(Tr);
%for i=1:n1
 %new_tr(i,:)=strrep(Tr(i,:),'X',' ');
%end
%lm=str2num(new_tr);marker=label(lm);chrgl=chr(lm);poslm=pos(lm);
%var_name_marker={'Bin';'Bin'
%Marker=marker;chromosome=chrgl;position=poslm;p_f=pr;beta=b(2:end);sebeta=SEb;r2=r2lz;
%opt_result=table(Marker,chromosome,position,f,p_f,beta,sebeta,r2);
%SNP=label;Chromosome=chr;Position=pos;Ft=ft;P_F=plm;Beta=effect;SEbeta=seblm;R2=r2lm;
%all_result=table(SNP,Chromosome,Position,Ft,P_F,Beta,SEbeta,R2);

if mdl==1
[n1,m1]=size(Tr);
for i=1:n1
 new_tr(i,:)=strrep(Tr(i,:),'X',' ');
end
lm=str2num(new_tr);marker=label(lm);chrgl=chr(lm);poslm=pos(lm);
Marker=marker;chromosome=chrgl;position=poslm;Ep_f=pr;Ebeta=b(2:end);Esebeta=SEb;Er2=r2lz;
opt_result=table(Marker,chromosome,position,f,Ep_f,Ebeta,Esebeta,Er2);
ESNP=label;Chromosome=chr;Position=pos;EFt=Eft;EP_F=Eplm;EBeta=Eeffect;ESEbeta=Eseblm;ER2=Er2lm;
all_result=table(ESNP,Chromosome,Position,EFt,EP_F,EBeta,ESEbeta,ER2);
lmbx=x(:,lm);
if gr==1
xlswrite(genotype,lmbx,1);
%xlswrite(genotype,Bx,2);
end
elseif mdl==2
%for i=1:p1
 %new_tr(i,:)=strsplit(TR(i,:),'S');
%end
%jm=1;jl=1;
% for im=jm:length(TR)
     % new_tr1(jl,:)=strsplit(TR(im,:),'S');
     %  jl=jl+1;
 %end
%lm_seq0=char(new_tr(:,2));lm_seq0=str2num(lm_seq0);lm_seq0(:,2)=lm_seq0;
%lm_seq1=char(new_tr1(:,2));lm_seq1=str2num(lm_seq1);
%lm_seq2=char(new_tr1(:,3));lm_seq2=str2num(lm_seq2);lm_seq3=[lm_seq1 lm_seq2];
%lm_seq(:,1:2)=lm_seq0;lm_seq(p1+1:length(TR),1:2)=lm_seq3;
%marker=label(lm_seq);chrgl=chr(lm_seq);startgl=[start(lm_seq) stop(lm_seq)];
%poslm=Binsize(lm_seq);
%clear lm_seq0 lm_seq1 lm_seq3 lm_seq2 new_str new_str1
%output the optmal result only the bin size
%if not the bin, then just only the pos, no bin size 
lm=trnum(tr,:);%Marker=lm;
   
%Marker=marker(lm);Chr=chrgl(lm);Position=poslm(lm);
Ep_f=pr;Ebeta=b(2:end);Esebeta=SEb;Er2=r2lz;
%opt_result=table(Marker(:,1),Marker(:,2),Chr(:,1),Chr(:,2),Position(:,1),Position(:,2),f,p_f,beta,sebeta,r2);
%opt_result.Properties.VariableNames={'Marker1' 'Marker2' 'Chr' 'Chr' 'Position' 'Position' 'f' 'Chr' 'p_f' 'beta' 'sebeta' 'r2'};
%opt_result=table(Marker,f,p_f,beta,sebeta,r2);
%save([opt_name '.mat'],'opt_result','-v7.3')
%output the all result
%SNP=marker;Chr=chrgl;Position=poslm;
EFt=Eft;EP_F=Eplm;EBeta=Eeffect;ESEbeta=Eseblm;ER2=Er2lm;ESNP=trnum;
%all_result=table(SNP,Ft,P_F,Beta,SEbeta,R2);
%save([allresult_name '.mat'],'all_result','-v7.3')
%clear marker chrgl poslm 
else
    error('Dosent chosen the right model')
end






