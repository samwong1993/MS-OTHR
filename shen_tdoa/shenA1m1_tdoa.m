function [Pim,ym,iter,optvalCurr] = shenA1m1_tdoa(ymPrev,tdoa_meas,x,M,N,c)
% PerfectP = 0: unknown permutation matrix
% metric = 1: l2-norm
iter=0;
optvalCurr=1;
optvalPrev=-1;
Omega = ones(N-1,N-1)+eye(N-1); inv_Omega =inv(Omega);
while abs(optvalCurr-optvalPrev)>1e-3
    cvx_begin quiet
        variables ym(2,M) Pim(M,M,N-1) mu1(N-1,M) tt(M);
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
                sum(Pim(m,:,i))==1;
                for j=1:M
                    Pim(m,j,i)<=1;
                    Pim(m,j,i)>=0;
                    sum(Pim(:,j,i))==1;        
                end
            end
        end
        for i = 1:M
            mu1(:,i)'*inv_Omega*mu1(:,i) <= tt(i)
        end
    cvx_end
    optvalPrev=optvalCurr;
    optvalCurr=cvx_optval;
    ymPrev=ym;
    iter=iter+1;
end

end