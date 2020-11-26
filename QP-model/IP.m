function [P_tau0 param] = IP(M,G,t,P_tau,K,param)
    G_t = G(1:M-1,1:M)*t;
    num = size(param.cut1,3);
    cut = '';
    for iter = 1:num
        for i = 1:M-2
            cut = cut + "trace(P"+string(i)+"'*param.cut"+string(i)+"(:,:,"+string(iter)+"))+";
        end
        i = i + 1;
        cut = cut + "trace(P"+string(i)+"'*param.cut"+string(i)+"(:,:,"+string(iter)+"))<=K*(M-1)-1;";
    end
    str = "";
    for i = 1:M-1
        str = str + "variable P"+string(i)+"("+string(K)+","+string(K)+") integer;";
    end
    cvx_begin quiet
    cvx_precision best
    eval(str)
    variable z(M-1,K)
    minimize sum(sum(z))
    for k = 1:M-1
        eval("P_P_tau"+string(k)+" = P_tau("+string(k)+",:)*P"+string(k)+";");
        for j = 1:K
            - z(k,j) <= G_t(k,j) - eval("P_P_tau"+string(k)+"(j)") <= z(k,j)
        end
    end
    for i = 1:K
        for j = 1:K
            for k = 1:M-1
                0 <= eval("P"+string(k)+"(i,j)") <= 1;
            end
        end
    end
    for i = 1:K
        for k = 1:M-1
            eval("sum(P"+string(k)+"(i,:)) == 1")
        end
    end
    for j = 1:K
        for k = 1:M-1
            eval("sum(P"+string(k)+"(:,j)) == 1")
        end
    end
    eval(cut);
    cvx_end
    P_tau0 = [];
    for i = 1:M-1
        P_tau0 = [P_tau0;eval("P_tau("+string(i)+",:)*P"+string(i))];
    end
    update = '';
    for i = 1:M-1
         update = update + "param.cut"+string(i)+"(:,:,num+1) = P"+string(i)+";";
    end
    eval(update);
    for i = 1:M-1
        eval("param.P"+string(i)+" = P"+string(i)+";");
    end
end





% G_t = G(1:M-1,1:M)*t;
% cvx_begin
% cvx_precision best
% variable P0(3,3) integer
% variable P1(3,3) integer
% variable P2(3,3) integer
% variable P3(3,3) integer
% variable z(M-1,3)
% minimize sum(sum(z))
% P_P_tau0 = P_tau(1,:)*P0;
% for j = 1:3
%     - z(1,j) <= G_t(1,j) - P_P_tau0(j) <= z(1,j)
% end
% P_P_tau1 = P_tau(2,:)*P1;
% for j = 1:3
%     - z(2,j) <= G_t(2,j) - P_P_tau1(j) <= z(2,j)
% end
% P_P_tau2 = P_tau(3,:)*P2;
% for j = 1:3
%     - z(3,j) <= G_t(3,j) - P_P_tau2(j) <= z(3,j)
% end
% P_P_tau3 = P_tau(4,:)*P3;
% for j = 1:3
%     - z(4,j) <= G_t(4,j) - P_P_tau3(j) <= z(4,j)
% end
% for i = 1:3
%     for j = 1:3
%         0 <= P0(i,j) <= 1;
%         0 <= P1(i,j) <= 1;
%         0 <= P2(i,j) <= 1;
%         0 <= P3(i,j) <= 1;
%     end
% end
% for i = 1:3
%     sum(P0(i,:)) == 1
%     sum(P1(i,:)) == 1
%     sum(P2(i,:)) == 1
%     sum(P3(i,:)) == 1
% end
% for j = 1:3
%     sum(P0(:,j)) == 1
%     sum(P1(:,j)) == 1
%     sum(P2(:,j)) == 1
%     sum(P3(:,j)) == 1
% end
% cvx_end
% P_tau0 = [P_tau(1,:)*P0;P_tau(2,:)*P1;P_tau(3,:)*P2;P_tau(4,:)*P3];