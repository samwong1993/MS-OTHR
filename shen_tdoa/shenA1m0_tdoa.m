function [Pim,ym,iter,optvalCurr] = shenA1m0_tdoa(ymPrev,tdoa_meas,x,M,N,c)
% metric = 0: universal upper bound

iter=0;
optvalCurr=1;
optvalPrev=-1;

while abs(optvalCurr-optvalPrev)>1e-3
    cvx_begin quiet
        variables ym(2,M) Pim(M,M,N-1) mu1;
        minimize mu1
        subject to
        for m=1:M
            ref = ymPrev(:,m)-x(:,1);
            refnorm = norm(ref);
            for i=1:N-1
                curr = ymPrev(:,m)-x(:,i+1);
                taylor = norm(curr) - refnorm + (curr/norm(curr) - ref/refnorm)'*(ym(:,m)-ymPrev(:,m));
                Pim(m,:,i)*tdoa_meas(i,:)'-(1/c)*taylor <= mu1;        
                Pim(m,:,i)*tdoa_meas(i,:)'-(1/c)*taylor >= -mu1;
                sum(Pim(m,:,i))==1;
                for j=1:M
                    Pim(m,j,i)<=1;
                    Pim(m,j,i)>=0;
                    sum(Pim(:,j,i))==1;        
                end
            end
        end
    cvx_end
    optvalPrev=optvalCurr;
    optvalCurr=cvx_optval;
    ymPrev=ym;
    iter=iter+1;
end

end