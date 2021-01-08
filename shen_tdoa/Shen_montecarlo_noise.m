%% Multiple sources
clc
clear all;
LOOP = 1; % Number of Monte Carlo simulation
MapConfig = 10; % "0": Large size; "1": Close sensors; else: Normal
InitialSel=1; % "0": generate a new initialization each iteration; "1": using the same initialization in each iteration
NoiseType = 0; % "0": Gaussian; else: Uniform
metric=1; % Algorithm 1: "0": universal upper bound; "1": l2-norm; "2": l1-norm
refineSel=0; % Way to determine the allocation; "0": Method in Paper; "1": Projecting to nearest permutation matrix
rmetric=1; % Algorithm 3: "0": universal upper bound; "1": l2-norm; "2": l1-norm
c=1;

if MapConfig == 0 % Large Size
    K=7; M=20;
    s=zeros(2,M); % anchor nodes
    s(:,1)=[-90;-90]; s(:,2)=[-90;-45]; s(:,3)=[-90;0]; s(:,4)=[-90;45]; s(:,5)=[-90;90];
    s(:,6:10)=s(:,1:5)+[[60;0] [60;0] [60;0] [60;0] [60;0]];
    s(:,11:15)=s(:,6:10)+[[60;0] [60;0] [60;0] [60;0] [60;0]];
    s(:,16:20)=s(:,11:15)+[[60;0] [60;0] [60;0] [60;0] [60;0]];
    xTrue=zeros(2,K); % source nodes
    xTrue(:,1)=[0;70]; xTrue(:,2)=[60;60]; xTrue(:,3)=[-60;30]; xTrue(:,4)=[-65;-40];
    xTrue(:,5)=[70;-60]; xTrue(:,6)=[0;-65]; xTrue(:,7)=[10;0];

elseif MapConfig==1 % Close Sources
    K=3; M=9;
    s=zeros(2,M); % anchor nodes
    s(:,1)=[40;40]; s(:,2)=[40;-40]; s(:,3)=[-40;40]; s(:,4)=[-40;-40]; s(:,5)=[40;0];
    s(:,6)=[0;40]; s(:,7)=[-40;0]; s(:,8)=[0;-40]; s(:,9)=[10;5];
    xTrue=zeros(2,K); % source nodes
    xTrue(:,1)=[10;-10]; xTrue(:,2)=[20;0]; xTrue(:,3)=[0;-10];
else
    
    %Attack example 1 
s = [40,40,-40,-40,40;40,-40,40,-40,0];
xTrue = [10,20,0;-10,0,-10];
%     %Attack example 2 
%     s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40];
%     xTrue = [10,-20,30;-100,-25,20];
%     %Attack example 3 
%     s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40];
%     xTrue = [10,-20,30;-10,-25,20];
% % 	%Attack example 4 
%     s = [800,300,1500,-600;600,450,-1000,700];
%     xTrue = [1000,-200,300;-1000,-250,200];
% %   example 4 
%     s = [-90,-90,-90,-90,-90,-30,-30,-30,-30,-30,30,30,30,30,30,90,90,90,90,90;-90,-45,0,45,90,-90,-45,0,45,90,-90,-45,0,45,90,-90,-45,0,45,90];
%     xTrue = [0,60,-60,-65,70,0,10;70,60,30,-40,-60,-65,0];

%     x=zeros(2,N); % anchor nodes
%     x(:,1)=[40;40]; x(:,2)=[40;-40]; x(:,3)=[-40;40]; x(:,4)=[-40;-40];
%     x(:,5)=[40;0]; x(:,6)=[0;40]; x(:,7)=[-40;0]; x(:,8)=[0;-40];
%     y=zeros(2,M); % source nodes
%     y(:,1)=[100;-100]; y(:,2)=[-200;-25]; y(:,3)=[30;200];
end
[d,M] = size(s); K = size(xTrue,2);
MFac=factorial(K); % for setting rmse vector
rmseInt=0; rmseFinal=0; rmsePer=0; % rmse of variables

varNos = [1 0.316227766016838 0.1 0.031622776601684 0.01 0.003162277660168 0.001];
SNR=10.*log10(1./varNos);
for idx_SNR = 3%1:length(SNR)
    for idx_seed = 1:70
    %% Generating measurements
%     rand('seed',idx_seed-1); randn('seed',idx_seed-1); % using the same set of random numbers
    %Generate noise
	G = zeros(M-1,M);
    G(:,1) = -ones(M-1,1);
    for i = 1:M-1
        G(i,i+1) = 1;
    end
    sigma_t = varNos(idx_SNR);
    noise_t0 = randn(M,K);
    noise_t = (sigma_t*G*noise_t0);
    tTrue=zeros(M,K); % true propagation time
    for i=1:M
        for j=1:K
            tTrue(i,j)=norm(s(:,i)-xTrue(:,j))/c;
        end
    end

    % Setting the first sensor as reference sensor, get the TDOA measurements
    delta_tTrue = zeros(M-1,K); tau = zeros(M-1,K);
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
%     sortIdx=zeros(K,M-1);
%     for i=1:M-1
%         [tau(i,:),sortIdx(:,i)]=sort(delta_tTrue(i,:),'ascend'); % Sort t by the ascending order of the TOA measurements
%     end
%     tau = tau + noise_t;
    %% Algorithm
    cvx_solver MOSEK
    if InitialSel==0
        ymPrev=randn(2,K);
    else
    %     if idx==1
            ymPrev=randn(2,K);
            ymPrevOne=ymPrev;
    %     else
    %        ymPrev=ymPrevOne; 
    %     end
    end

    %% Algorithm 1
    if metric==0
        [Pim,ym,iter,optvalCurr] = shenA1m0_tdoa(ymPrev,tau,s,K,M,c);
    elseif metric==1
        [Pim,ym,iter,optvalCurr] = shenA1m1_tdoa(ymPrev,tau,s,K,M,c);
    elseif metric==2
        [Pim,ym,iter,optvalCurr] = shenA1m2_tdoa(ymPrev,tau,s,K,M,c);
    end
    yE=ym; % Source estimate
    rmse_new = rmse(yE,xTrue,MFac,K);
    rmseInt = rmseInt + rmse_new;
    yIni = ym;

    %% Algorithm 2
    if refineSel==0
        for i=1:M-1
            Z=ones(K,K);
            for im=1:K
            [maxval1,maxidx1]=max(Pim(:,:,i).*Z);
            [maxval2,maxidx2]=max(maxval1);
            Pim(:,maxidx2,i)=zeros(K,1);
            Pim(maxidx1(maxidx2),:,i)=zeros(1,K);
            Pim(maxidx1(maxidx2),maxidx2,i)=1;
            Z(:,maxidx2)=zeros(K,1);
            Z(maxidx1(maxidx2),:)=zeros(1,K);
            end
        end
    else
        IndexAll=perms([1:K]);
        % Transforming the ordering of IndexAll to permutation matrices
        PiAll=zeros(K,K,MFac);
        for m=1:MFac
            for m1=1:K
                PiAll(m1,IndexAll(m,m1),m)=1;
            end
        end
        % Projecting the matrix Pim to the nearest permutation matrix
        mseP=zeros(MFac,M-1);
        for i=1:M-1
            mseMin=1e4;
            for m=1:MFac
                mseP(m,i)=norm(PiAll(:,:,m)-Pim(:,:,i),'fro')^2;        
            end
        end
        [mseSort,IdxSort]=sort(mseP);
        Pim1=zeros(K,K,M-1); % refined permutation matrix
        for i=1:M-1
            Pim1(:,:,i)=PiAll(:,:,IdxSort(1,i));
        end
        Pim=Pim1;
    end

    %% Algorithm 3
    if rmetric==0
        [yE,iterR,optvalCurr] = shenA3m0_tdoa(ym,Pim,tau,s,K,M,c);
    elseif rmetric==1
        [yE,iterR,optvalCurr] = shenA3m1_tdoa(ym,Pim,tau,s,K,M,c);
    elseif rmetric==2
        [yE,iterR,optvalCurr] = shenA3m2_tdoa(ym,Pim,tau,s,K,M,c);
    end

    rmseRe = rmse(yE,xTrue,MFac,K);    
    rmseFinal = rmseFinal + rmseRe;
    rmseCurrent=sqrt(rmseFinal/idx_seed/K);
    [rmseCurrent idx_seed];

    [P] = compute_err(yE',xTrue);
    x_s = P*yE';
    for i = 1:size(xTrue,2)
        err(i) = norm(x_s(i,:) - xTrue(:,i)');
    end
    err
    fid=fopen("shen_model_1_M5_SNR"+string(SNR(idx_SNR))+".txt","a+");
    fprintf(fid,"%2.4f",err(1));
    for i = 2:size(xTrue,2)
        fprintf(fid,",%2.4f",err(i));
    end
    fprintf(fid,"\n");
    fclose(fid);
    end
end
% % Print out result
% rmseInt=sqrt(rmseInt/LOOP/M);
% rmseFinal=sqrt(rmseFinal/LOOP/M);
% rmsePer=sqrt(rmsePer/LOOP/M);
% iterInt=iter/LOOP;
% iterFinal=iterR/LOOP;
% mean(iterFinal);
% plot(x(1,:), x(2,:),'bo', y(1,:), y(2,:), 'r*',yE(1,:), yE(2,:), 'ks');
% legend({'Sensors', 'Sources', 'Recovered Sources'});