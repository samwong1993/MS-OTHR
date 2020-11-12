function [P] = compute_err(x_rec,emitter)
    K = size(emitter,2);
    cvx_begin quiet
    cvx_precision best
    variable z(K,3)
    variable P(K,K) integer
    minimize sum(sum(z))
    for i = 1:K
        for j = 1:3
            - z(i,j) <= P*x_rec - emitter' <= z(i,j)
        end
    end
    for i = 1:K
        for j = 1:K
            0 <= P(i,j) <= 1
        end
    end
    cvx_end
end



