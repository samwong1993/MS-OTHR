varNos = [1 0.316227766016838 0.1 0.031622776601684 0.01 0.003162277660168 0.001];
SNR=10.*log10(1./varNos);
err = [];
clf
for idx_SNR = 1:length(SNR)
    filename = ".\model_1_SNR"+string(SNR(idx_SNR))+".txt";
    [err_1,err_2,err_3]=textread(filename,'%f%f%f','delimiter',',');
    err = [err (mean(err_1) + mean(err_2) + mean(err_3))];
end
err_shen = [];
for idx_SNR = 1:length(SNR)
    filename = ".\shen_model_1_SNR"+string(SNR(idx_SNR))+".txt";
    [err_1,err_2,err_3]=textread(filename,'%f%f%f','delimiter',',');
    err_shen = [err_shen (mean(err_1) + mean(err_2) + mean(err_3))];
end

plot(SNR,err_shen,'^-', 'linewidth', 1.1, 'markerfacecolor', [29, 191, 151]/255);
hold on
plot(SNR,err,'ok-','linewidth',1.1,'markerfacecolor',[36, 169, 225]/255);

% for i = 1:m
%     plot(sigma,10^3*error(i,:),'b--o') 
% end
xlabel('SNR (dB)')
ylabel('RMSE (m)')
set(gca,'yscale','log')

%Attack example 1 
s = [40,40,-40,-40,40,0,-40,0,10;40,-40,40,-40,0,40,0,-40,5];
xTrue = [10,20,0;-10,0,-10];
% % %Attack example 2 
% s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40];
% xTrue = [10,-20,30;-100,-25,20];
% %Attack example 3 
% s = [40,40,-40,-40,40,0,-40,0;40,-40,40,-40,0,40,0,-40];
% xTrue = [10,-20,30;-10,-25,20];
% %Attack example 4 
% s = [800,300,1500,-600;600,450,-1000,700];
% xTrue = [1000,-200,300;-1000,-250,200];

[d,M] = size(s); K = size(xTrue,2);
Omega = ones(M-1,M-1)+eye(M-1);
sigma = [1 0.316227766016838 0.1 0.031622776601684 0.01 0.003162277660168 0.001];
crlb = zeros(length(sigma),1);
index = [1:M] ;
for idx = 1:K
    for noise_level = 1:length(sigma)
        eval("crlb(noise_level) = crlb(noise_level) + TDOALocCRLB(s,xTrue(:,"+string(idx)+"),Omega,sigma(noise_level));");
    end
end
SNR=10.*log10(1./varNos);
plot(SNR,crlb,'r','linewidth',1.5);
legend("Shen's Alg with L-2 norm",'AMLC','CRLB');