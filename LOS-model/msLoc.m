function [t_sum,obj_sum,x] = msLoc(s,P_tau0,Omega,inv_Omega,d,K,M,c)
% Input:
% s: (d-by-M) sensor locations
% P_tau: (M-1,K) Permuted time-difference
% Omega: (M-1)-by-(M-1) covariance matrix
% d: dimension; K: number of sources; M: number of anchors,
% c: speed of light

% Output:
% t_sum: (M-by-K) time-of-flight
% obj_sum: (K-by-1) objective value of k-th source

t_sum = zeros(M,K); obj_sum = zeros(K,1); x = zeros(d,K);

G = zeros(M-1,M); G(:,1) = -ones(M-1,1);
for i = 1:M-1
    G(i,i+1) = 1;
end

for k = 1:K
    x(:,k) = TDOALoc(s, c*P_tau0(:,k), Omega);
    for i = 1:M
       t_sum(i,k) = norm(x(:,k) - s(:,i))/c;
    end
    obj_sum(k) = (G*t_sum(:,k) - P_tau0(:,k))'*inv_Omega*(G*t_sum(:,k) - P_tau0(:,k));

end