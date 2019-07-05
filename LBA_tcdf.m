function F = LBA_tcdf(t, A, b, v, sv)
% Get CDF of first passage time of ith accumulator in LBA model
% F = LBA_tcdf(t, A, b, v, sv)


g = (b-A-t.*v)./(t.*sv); % chizumax
h = (b-t.*v)./(t.*sv);  % chizu
i = b-t.*v; % chiminuszu
j = i-A; % xx

g=real(g);
h=real(h);
i=real(i);
j=real(j);

temp1=real(normcdf(g));
temp2=real(normcdf(h));

id=temp1>0.9999;
temp1(id,1)=0.9999;
id=temp1<=0.0001;
temp1(id,1)=0.0001;

id=temp2>0.9999;
temp2(id,1)=0.9999;
id=temp2<0.0001;
temp2(id,1)=0.0001;

tmp1 = t.*sv.*(normpdf(g)-normpdf(h));
tmp2 = j.*temp1 - i.*temp2;

F = 1 + (tmp1 + tmp2)./A;