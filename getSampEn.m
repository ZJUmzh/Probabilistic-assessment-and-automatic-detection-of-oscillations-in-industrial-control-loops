function [sampEn,phi] = getSampEn(x,varargin)
narginchk(1,9)
x = squeeze(x);
x = x(:);
p = inputParser;
chk1 = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
chk2 = @(x) isscalar(x) && (x > 0);
addRequired(p,'x',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'m',2,chk1);
addParameter(p,'tau',1,chk1);
addParameter(p,'r',.2*std(x,1),chk2);
addParameter(p,'logx',exp(1),chk2);
parse(p,x,varargin{:})
m = p.Results.m; tau = p.Results.tau; 
r = p.Results.r; logx = p.Results.logx; 
method = 'single'; 
switch method
    case 'single'
        x = single(x);
        r = single(r);
end
r = r*std(x);
lenx = length(x);
phi = zeros(1,2,method);
for p = 1:2
    rows = m+p-1;
    cols = lenx-(rows-1)*tau;
    dataMat = zeros(rows,cols,method);
    if rows < cols
        for ii = 1:rows
            dataMat(ii,:) = x((ii-1)*tau+1:(ii-1)*tau+cols);
        end
    else
        for ii = 1:cols
            dataMat(:,ii) = x(ii:tau:ii+rows*tau-1);
        end
    end  
    dataMat = dataMat(:,1:lenx-m*tau);
    d = pdist(dataMat','chebychev');
    phi(p) = sum(d<=r);
end
sampEn = log(phi(1)/phi(2))/log(logx);

end