%% Multiple sources
clc
clear all;
LOOP = 1; % Number of Monte Carlo simulation
MapConfig = 10; % "0": Large size; "1": Close sensors; else: Normal
varNos = 1; % variance of noise
SNR=10*log10(1/varNos);
InitialSel=1; % "0": generate a new initialization each iteration; "1": using the same initialization in each iteration
NoiseType = 0; % "0": Gaussian; else: Uniform
metric=1; % Algorithm 1: "0": universal upper bound; "1": l2-norm; "2": l1-norm
refineSel=0; % Way to determine the allocation; "0": Method in Paper; "1": Projecting to nearest permutation matrix
rmetric=0; % Algorithm 3: "0": universal upper bound; "1": l2-norm; "2": l1-norm
c=1;

if MapConfig == 0 % Large Size
    M=7; N=20;
    x=zeros(2,N); % anchor nodes
    x(:,1)=[-90;-90]; x(:,2)=[-90;-45]; x(:,3)=[-90;0]; x(:,4)=[-90;45]; x(:,5)=[-90;90];
    x(:,6:10)=x(:,1:5)+[[60;0] [60;0] [60;0] [60;0] [60;0]];
    x(:,11:15)=x(:,6:10)+[[60;0] [60;0] [60;0] [60;0] [60;0]];
    x(:,16:20)=x(:,11:15)+[[60;0] [60;0] [60;0] [60;0] [60;0]];
    y=zeros(2,M); % source nodes
    y(:,1)=[0;70]; y(:,2)=[60;60]; y(:,3)=[-60;30]; y(:,4)=[-65;-40];
    y(:,5)=[70;-60]; y(:,6)=[0;-65]; y(:,7)=[10;0];

elseif MapConfig==1 % Close Sources
    M=3; N=9;
    x=zeros(2,N); % anchor nodes
    x(:,1)=[40;40]; x(:,2)=[40;-40]; x(:,3)=[-40;40]; x(:,4)=[-40;-40]; x(:,5)=[40;0];
    x(:,6)=[0;40]; x(:,7)=[-40;0]; x(:,8)=[0;-40]; x(:,9)=[10;5];
    y=zeros(2,M); % source nodes
    y(:,1)=[10;-10]; y(:,2)=[20;0]; y(:,3)=[0;-10];
else
    M=3; N=8; % Normal Sources
    
%     %Attack example 1 
%     s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40];
%     xTrue = [100,-200,30;-100,-25,200];
    %Attack example 2 
    s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40];
    xTrue = [10,-20,30;-100,-25,20];
%     %Attack example 3 
%     s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40];
%     xTrue = [10,-20,30;-10,-25,20];
% 	%Attack example 4 
%     s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40];
%     xTrue = [100,-200,-100,100];

    x = s;
    y = xTrue;
    
%     x=zeros(2,N); % anchor nodes
%     x(:,1)=[40;40]; x(:,2)=[40;-40]; x(:,3)=[-40;40]; x(:,4)=[-40;-40];
%     x(:,5)=[40;0]; x(:,6)=[0;40]; x(:,7)=[-40;0]; x(:,8)=[0;-40];
%     y=zeros(2,M); % source nodes
%     y(:,1)=[100;-100]; y(:,2)=[-200;-25]; y(:,3)=[30;200];
end

MFac=factorial(M); % for setting rmse vector
rmseInt=0; rmseFinal=0; rmsePer=0; % rmse of variables

varNos = [1 0.316227766016838 0.1 0.031622776601684 0.01 0.003162277660168 0.001];
SNR=10.*log10(1./varNos);
for idx_SNR = 1:length(SNR)
    for idx_seed = 1:30
    %% Generating measurements
    rand('seed',idx_seed-1); randn('seed',idx_seed-1); % using the same set of random numbers
    if NoiseType==0 % Gaussian noise
        nm=randn(N,M)*sqrt(varNos(idx_SNR));
    else % uniform noise
        b=sqrt(3*varNos(idx_SNR));
        a=-b;
        nm=a + (b-a).*rand(N,M);
    end
    t0=zeros(N,M); % true propagation time
    for i=1:N
        for j=1:M
            t0(i,j)=norm(x(:,i)-y(:,j))/c+nm(i,j);
        end
    end
    % Setting the first sensor as reference sensor, get the TDOA measurements
    tdoa = zeros(N-1,M); tdoa_meas = zeros(N-1,M);
    for i = 1:N-1
        tdoa(i,:) = t0(i+1,:) - t0(1,:);
    end
    sortIdx=zeros(M,N-1);
    for i=1:N-1
        [tdoa_meas(i,:),sortIdx(:,i)]=sort(tdoa(i,:),'ascend'); % Sort t by the ascending order of the TOA measurements
    end

    %% Algorithm
    cvx_solver gurobi
    if InitialSel==0
        ymPrev=randn(2,M);
    else
    %     if idx==1
            ymPrev=randn(2,M);
            ymPrevOne=ymPrev;
    %     else
    %        ymPrev=ymPrevOne; 
    %     end
    end

    %% Algorithm 1
    if metric==0
        [Pim,ym,iter,optvalCurr] = shenA1m0_tdoa(ymPrev,tdoa_meas,x,M,N,c);
    elseif metric==1
        [Pim,ym,iter,optvalCurr] = shenA1m1_tdoa(ymPrev,tdoa_meas,x,M,N,c);
    elseif metric==2
        [Pim,ym,iter,optvalCurr] = shenA1m2_tdoa(ymPrev,tdoa_meas,x,M,N,c);
    end
    yE=ym; % Source estimate
    rmse_new = rmse(yE,y,MFac,M);
    rmseInt = rmseInt + rmse_new;
    yIni = ym;

    %% Algorithm 2
    if refineSel==0
        for i=1:N-1
            Z=ones(M,M);
            for im=1:M
            [maxval1,maxidx1]=max(Pim(:,:,i).*Z);
            [maxval2,maxidx2]=max(maxval1);
            Pim(:,maxidx2,i)=zeros(M,1);
            Pim(maxidx1(maxidx2),:,i)=zeros(1,M);
            Pim(maxidx1(maxidx2),maxidx2,i)=1;
            Z(:,maxidx2)=zeros(M,1);
            Z(maxidx1(maxidx2),:)=zeros(1,M);
            end
        end
    else
        IndexAll=perms([1:M]);
        % Transforming the ordering of IndexAll to permutation matrices
        PiAll=zeros(M,M,MFac);
        for m=1:MFac
            for m1=1:M
                PiAll(m1,IndexAll(m,m1),m)=1;
            end
        end
        % Projecting the matrix Pim to the nearest permutation matrix
        mseP=zeros(MFac,N-1);
        for i=1:N-1
            mseMin=1e4;
            for m=1:MFac
                mseP(m,i)=norm(PiAll(:,:,m)-Pim(:,:,i),'fro')^2;        
            end
        end
        [mseSort,IdxSort]=sort(mseP);
        Pim1=zeros(M,M,N-1); % refined permutation matrix
        for i=1:N-1
            Pim1(:,:,i)=PiAll(:,:,IdxSort(1,i));
        end
        Pim=Pim1;
    end

    %% Algorithm 3
    if rmetric==0
        [yE,iterR,optvalCurr] = shenA3m0_tdoa(ym,Pim,tdoa_meas,x,M,N,c);
    elseif rmetric==1
        [yE,iterR,optvalCurr] = shenA3m1_tdoa(ym,Pim,tdoa_meas,x,M,N,c);
    elseif rmetric==2
        [yE,iterR,optvalCurr] = shenA3m2_tdoa(ym,Pim,tdoa_meas,x,M,N,c);
    end

    rmseRe = rmse(yE,y,MFac,M);    
    rmseFinal = rmseFinal + rmseRe;
    rmseCurrent=sqrt(rmseFinal/idx_seed/M);
    [rmseCurrent idx_seed];

    
    cvx_solver gurobi
    [P] = compute_err(yE',y);
    x_s = P*yE';
    for i = 1:size(y,2)
        err(i) = norm(x_s(i,:) - y(:,i)');
    end
    err
    fid=fopen("shen_model_2_SNR"+string(SNR(idx_SNR))+".txt","a+");
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