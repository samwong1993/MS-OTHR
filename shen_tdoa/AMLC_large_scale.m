%% initialization
% addpath("C:\Users\ympun\Downloads\cvx-w64")
% cvx_setup

clear all
clc;
c = 1;
% range = 100;initial = -50;
% %dim:d=2 number of sources:K=2 senors M=8;
% d = 2; K = 5 ; M = 12; eta = 0.01;
% xTrue = initial + range*rand(d,K);
% s = initial + range*rand(d,M);   % Sensor Location
% rng(23)
% d = 2; K = 7; M = 20; 
% range = 1000; initial = -500;
% xTrue = initial + range*rand(d,K);
% s = initial + range*rand(d,M);   % Sensor Location
% Example 1 (Shen)
% s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40]; xTrue = [10,-20,30;-10,-25,20];
% % Example 6 (Shen) --- Large-scale
% s = [-90,-90,-90,-90,-90,-30,-30,-30,-30,-30,30,30,30,30,30,90,90,90,90,90;-90,-45,0,45,90,-90,-45,0,45,90,-90,-45,0,45,90,-90,-45,0,45,90];
% xTrue = [0,60,-60,-65,70,0,10;70,60,30,-40,-60,-65,0];
%Attack example 1 
s = [40,40,-40,-40,40,0,-40,0,10;40,-40,40,-40,0,40,0,-40,5];
xTrue = [10,20,0;-10,0,-10];
%Attack example 2 
% s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40];
% xTrue = [10,-20,30;-100,-25,20];
% %Attack example 3 
% s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40];
% xTrue = [10,-20,30;-10,-25,20];
%Attack example 4 
% s = [800,300,1500,-600;600,450,-1000,700];
% xTrue = [1000,-200,300;-1000,-250,200];
% %Attack example 5
% s = [800,300,1500,-600,1500,1000,-1000,-1000,0;600,450,-1000,700,800,0,0,-1000,-1500];
% xTrue = [1000,-200,300;-1000,-250,200];
%Attack example 6 
s = [0,0,0,0,0,0,0,0,0,0,0,-100,100,-200,200,-300,300,-400,400,-500,500;
    0,-100,100,-200,200,-300,300,-400,400,-500,500,0,0,0,0,0,0,0,0,0,0];
xTrue = [1500,500;2000,1000];
%Attack example 7 
s = [0,500,500,-500,-500;
    0,500,-500,500,-500];
xTrue = [3000,1000,6000,0,4000,5000,-3000,-1000;2000,5000,0,8000,6000,-2000,3000,-1000];
% xTrue(:,1) = [];
% Example 7 (Shen)
% s = [40,40,-40,-40,40,0,-40,0,10;40,-40,40,-40,0,40,0,-40,0]; xTrue = [10,20,0;-10,0,-10];
[d,M] = size(s); K = size(xTrue,2);
Omega = ones(M-1,M-1)+eye(M-1); inv_Omega =inv(Omega); % covariance matrix


% figure(1)
% hold on
varNos = [1 0.316227766016838 0.1 0.031622776601684 0.01 0.003162277660168 0.001];
SNR=10.*log10(1./varNos);
for idx_SNR = 7%1:length(SNR)
        
    for idx_seed = 1:5
    %% Generating measurements
%     rand('seed',idx_seed-1); randn('seed',idx_seed-1); % using the same set of random numbers
	%Generate noise
    clear param
    K = size(xTrue,2);
	G = zeros(M-1,M);
    G(:,1) = -ones(M-1,1);
    for i = 1:M-1
        G(i,i+1) = 1;
    end
    sigma_t = varNos(idx_SNR);
    noise_t0 = randn(M,K);
    noise_t = (sigma_t*G*noise_t0);
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
    for i=1:M-1
        orderPerm(:,i) = randperm(K)';
        PTrue(:,:,i) = I(orderPerm(:,i),:);
        tau(i,:) = delta_tTrue(i,:)* PTrue(:,:,i) + noise_t(i,:);
    end
	P_tau = tau; t = zeros(M,K); 
    obj_best = 9999999;
    %% solve IP
    ini = '';
    for i = 1:M-1
         ini = ini + "param.cut"+string(i)+"(:,:,1) = zeros(K,K);";
    end
    eval(ini);
    x_rec = [];
    for iter = 1:200
        if K ~= 1
            [P_tau0, param] = IP_los(G,param,K,M,t,P_tau);
        end
        obj = trace((G(1:M-1,1:M)*t - P_tau0)'*inv_Omega*(G(1:M-1,1:M)*t - P_tau0));
        fprintf("iter:%d obj:%2.8f K:%d\n",iter,obj,K);
        [t_sum,obj_sum,location] = msLoc(s,P_tau0,Omega,inv_Omega,d,K,M,c);
        t = t_sum;
        index = find(obj_sum >= 1e-1);
        dif_index = find(obj_sum <= 1e-1);
        if isempty(index)
            %Record the location of the last iterartion and break
            for i = 1:length(dif_index)
                eval("x_rec = [x_rec,location(:,"+string(dif_index(i))+")];");
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
                eval("x_rec = [x_rec,location(:,"+string(dif_index(i))+")];");
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
    fid=fopen("model_7_M5_SNR"+string(SNR(idx_SNR))+".txt","a+");
    fprintf(fid,"%2.4f",err(1));
    for i = 2:size(xTrue,2)
        fprintf(fid,",%2.4f",err(i));
    end
    fprintf(fid,"\n");
    fclose(fid);
%     plot(s(1,:), s(2,:),'bo', xTrue(1,:), xTrue(2,:), 'r*',x_rec(1,:), x_rec(2,:), 'ks');
%     legend({'Sensors', 'Sources', 'Recovered Sources'});
    end
end