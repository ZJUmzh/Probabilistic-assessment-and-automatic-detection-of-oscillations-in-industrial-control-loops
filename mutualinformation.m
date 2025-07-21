function MI=mutualinformation(u1,u2)
wind_size=floor(power(length(u1),1/3)+0.5);
x=[u1,u2];
n=wind_size;
[xrow,xcol]=size(x);
bin=zeros(xrow,xcol);
pmf=zeros(n,2);
for i=1:2
    minx=min(x(:,i));
    maxx=max(x(:,i));
    binwidth=(maxx-minx)/n;
    edges=minx+binwidth*(0:n);
    histcedges=[-inf edges(2:end-1) inf];
    [occur,~,bin(:,i)]=histcounts(x(:,i),histcedges);
    pmf(:,i)=occur(1:n)./xrow;
end
jointoccur=accumarray(bin,1,[n,n]);
jointpmf=jointoccur./xrow;
Hx=-(pmf(:,1))'*log2(pmf(:,1)+eps);
Hy=-(pmf(:,2))'*log2(pmf(:,2)+eps);
Hxy=-(jointpmf(:))'*log2(jointpmf(:)+eps);
MI=Hx+Hy-Hxy;
end