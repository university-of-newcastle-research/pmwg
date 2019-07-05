function pdf = LBA_n1PDF_reparam_real(t, log_A, log_b, log_v, sv,log_tau)
% Generates defective PDF for responses on node #1 (ie. normalised by
% probability of this node winning race)
%
% pdf = LBA_n1PDF(t, A, b, v, sv)
%


N = size(log_v,2);
A = exp(log_A);
b = exp(log_b);
v = exp(log_v);
tau = exp(log_tau);
if N > 2
    for i = 2:N
        tmp(:,i-1) = LBA_tcdf(t-tau,A,b,v(:,i),sv);
    end
    G = prod(1-tmp,2);
else
    G = 1-LBA_tcdf(t-tau,A,b,v(:,2),sv);
end
%pdf = G.*LBA_tpdf(t-tau,A,b,v(:,1),sv);
id=G<=10^-50;
G(id,1)=10^-50;

H=LBA_tpdf(t-tau,A,b,v(:,1),sv);
id=H<=10^-50;
H(id,1)=10^-50;

pdf=exp(log(G)+log(H));