function [R,C] =topsis(result, method)
%load allresult.mat
% [r,c]=topsis('allresult.mat','high')
% [r,c]=topsis('result.mat','high')
format bank
load(result)
[n,m]=size(x);
%original data needs to be converged
if method == "low"
    x=1./(x);
elseif  method == "neutral"
        for i=1:m
            X(:,i)=median(x(:,i))/(median(x(:,i))+abs(x(:,i)-median(x(:,i))));
        end
        x=X;
        [n,m]=size(x);
end
zh=zeros(1,m);
d1=zeros(1,n); %min matrix 
d2=zeros(1,n); %max matrix 
c=zeros(1,n);
%normalize
for i=1:m
    for j=1:n
        zh(i)=zh(i)+x(j,i)^2;
    end
end
for i=1:m
    for j=1:n
       x(j,i)=x(j,i)/sqrt( zh(i));
    end
end
%calculate distance 
xx=min(x);
dd=max(x);
for i=1:n
    for j=1:m
        d1(i)=d1(i)+(x(i,j)-xx(j))^2;
    end
    d1(i)=sqrt(d1(i));
end
for i=1:n
    for j=1:m
        d2(i)=d2(i)+(x(i,j)-dd(j))^2;
    end
    d2(i)=sqrt(d2(i));
end
%Calculate the closeness 
for i=1:n
    c(i)=d1(i)/(d2(i)+d1(i));
end
c'
[a,b]=sort(c,'descend');
C=a';
R=b';
end
