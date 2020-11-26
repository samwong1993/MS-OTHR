%% initialization
% addpath("C:\Users\ympun\Downloads\cvx-w64")
% cvx_setup

clear all
clc;
c = 10;
% rng(22);
range = 100;initial = -50;
%dim:d=2 number of sources:K=2 senors M=8;
d = 2; K = 3; M = 12; eta = 0.01;
xTrue = initial + range*rand(d,K);
s = initial + range*rand(d,M);   % Sensor Location
Omega = ones(M-1,M-1)+eye(M-1); inv_Omega =inv(Omega); % covariance matrix

tTrue = zeros(M,K);
for i = 1:M
    for k = 1:K
        tTrue(i,k) = norm(s(:,i)-xTrue(:,k)) / c;
    end
end

delta_tTrue = zeros(M-1,K);
for i = 1:M-1
    delta_tTrue(i,:) = tTrue(i+1,:) - tTrue(1,:);
end

PTrue = zeros(K,K,M-1); I = eye(K);
tau = zeros(M-1,K); orderPerm = zeros(K, M-1);
% t_noise=eta/2*randn(K,M);
t_noise = eta/2*rand(M,K);
 for i=1:M-1
     orderPerm(:,i) = randperm(K)';
     PTrue(:,:,i) = I(orderPerm(:,i),:);
     t_error = t_noise(i+1,:) - t_noise(1,:);
     tau(i,:) = (delta_tTrue(i,:) + t_error)* PTrue(:,:,i);
 end
 
 G = zeros(M-1,M);
 G(:,1) = -ones(M-1,1);
for i = 1:M-1
    G(i,i+1) = 1;
end

%% solve £¨IP£©
ini = '';
for i = 1:M-1
     ini = ini + "param.cut"+string(i)+"(:,:,1) = zeros(K,K);";
end
eval(ini);

P_tau = tau; t = zeros(M,K); x_rec = [];
for iter = 1:100
    if K ~= 1
        [P_tau0, param] = IP_los(G,param,K,M,t,P_tau);
    end
    obj = trace((G*t - P_tau0)'*inv_Omega*(G*t - P_tau0));
    fprintf("obj:%2.8f K:%d\n",obj,K);
    [t_sum,obj_sum,location] = msLoc(s,P_tau0,Omega,inv_Omega,d,K,M,c);
    t = t_sum;
    index = find(obj_sum >= 1e-2);
    dif_index = find(obj_sum <= 1e-2);
    if isempty(index)
        %Record the location of the last iterartion and break
        for i = 1:length(dif_index)
            x_rec = [x_rec,location(:,dif_index(i))];
        end
        break
    end
    if ~isempty(dif_index)
        %Update P_tau and K
        for i = 1:length(dif_index)
            P_tau = del(P_tau,P_tau0(:,dif_index(i)));
        end
        K = size(P_tau,2);
        %Record the location
        for i = 1:length(dif_index)
            x_rec = [x_rec,location(:,dif_index(i))];
        end
        t = t_sum(:,index);
        clear param
        ini = '';
        for i = 1:M-1
             ini = ini + "param.cut"+string(i)+"(:,:,1) = zeros(K,K);";
        end
        eval(ini);
    end
end
%Compute the localization error
[P] = compute_err(x_rec',xTrue);
x = P*x_rec';
for i = 1:size(xTrue,2)
    err(i) = norm(x(i,:) - xTrue(:,i)');
end
%Save the results
fid=fopen("results.txt","a+");
fprintf(fid,"%2.4f",err(1));
for i = 2:size(xTrue,2)
    fprintf(fid,",%2.4f",err(i));
end
fprintf(fid,"\n");
fclose(fid);