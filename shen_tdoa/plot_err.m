varNos = [1 0.316227766016838 0.1 0.031622776601684 0.01 0.003162277660168 0.001];
SNR=10.*log10(1./varNos);
err = [];
for idx_SNR = 1:length(SNR)
    filename = ".\model_3_SNR"+string(SNR(idx_SNR))+".txt";
    [err_1,err_2,err_3]=textread(filename,'%f%f%f','delimiter',',');
    err = [err (mean(err_1) + mean(err_2) + mean(err_3))/3];
end
err_shen = [];
for idx_SNR = 1:length(SNR)
    filename = ".\shen_model_3_SNR"+string(SNR(idx_SNR))+".txt";
    [err_1,err_2,err_3]=textread(filename,'%f%f%f','delimiter',',');
    err_shen = [err_shen (mean(err_1) + mean(err_2) + mean(err_3))/3];
end

plot(SNR,err_shen,'^-', 'linewidth', 1.1, 'markerfacecolor', [29, 191, 151]/255);
hold on
plot(SNR,err,'ok-','linewidth',1.1,'markerfacecolor',[36, 169, 225]/255);
legend("Shen's Alg with L-2 norm",'AMLC');
% for i = 1:m
%     plot(sigma,10^3*error(i,:),'b--o') 
% end
xlabel('SNR (dB)')
ylabel('RMSE (m)')
