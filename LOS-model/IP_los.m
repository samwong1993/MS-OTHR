function [P_tau, param] = IP_los(G,param,K,M,t,tau)
% constraints about the cutting set
e = ones(K,1);
num = size(param.cut1,3);
cut = '';
for iter = 1:num
    for i = 1:M-2
        cut = cut + "trace(P(:,:,"+string(i)+")'*param.cut"+string(i)+"(:,:,"+string(iter)+"))+";
    end
    i = i + 1;
    cut = cut + "trace(P(:,:,"+string(i)+")'*param.cut"+string(i)+"(:,:,"+string(iter)+"))<=K*(M-1)-1;";
end

% cvx_setup;
cvx_solver gurobi_2;
cvx_begin quiet
cvx_precision best
variable P(K,K,M-1) binary
variable z(M-1,K)

minimize sum(sum(z))
subject to

for i=1:M-1
    for k=1:K
        -z(i,k) <= G(i,:)*t(:,k) - tau(i,:)*P(k,:,i)' <= z(i,k);
    end
    P(:,:,i)'* e == e;
    P(:,:,i) * e == e;
end
eval(cut);
cvx_end

P_tau = zeros(M-1,K);
for i = 1:M-1
    P_tau(i,:) = tau(i,:)*P(:,:,i)';
end
% update the cutting set
update = '';
for i = 1:M-1
    update = update + "param.cut"+string(i)+"(:,:,num+1) = P(:,:,"+string(i)+");";
end
eval(update);
param.z = z;
end
