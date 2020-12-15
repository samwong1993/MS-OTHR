function [ym,iterR,optvalCurr] = shenA3m1_tdoa(ym,Pim,tdoa_meas,x,M,N,c)
% metric = 1: l2-norm
iterR=zeros(M,1);
ymPrev=ym;      
optvalCurr=1;
optvalPrev=-1;
Omega = ones(N-1,N-1)+eye(N-1); inv_Omega =inv(Omega);
while abs(optvalCurr-optvalPrev)>1e-3
	cvx_begin quiet
        variables ym(2,M) mu1(N-1,M) tt(M);
        minimize sum(tt)
        subject to
        for m=1:M
            ref = ymPrev(:,m)-x(:,1);
            refnorm = norm(ref);
            for i=1:N-1
                curr = ymPrev(:,m)-x(:,i+1);
                taylor = norm(curr) - refnorm + (curr/norm(curr) - ref/refnorm)'*(ym(:,m)-ymPrev(:,m));
                Pim(m,:,i)*tdoa_meas(i,:)'-(1/c)*taylor <= mu1(i,m);        
                Pim(m,:,i)*tdoa_meas(i,:)'-(1/c)*taylor >= -mu1(i,m);
            end
        end
        for i = 1:M
            mu1(:,i)'*inv_Omega*mu1(:,i)<=tt(i)
        end
    cvx_end
    optvalPrev=optvalCurr;
    optvalCurr=cvx_optval;
    ymPrev=ym;
    iterR(m)=iterR(m)+1;
end

end