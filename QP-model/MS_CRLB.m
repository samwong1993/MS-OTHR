sigma = [0 100 200 300 400 500 600 700 800 900 1000];
crlb = zeros(length(sigma),1);
M = 6;
K = 4;
index = [1 2 3 4 5 6] ;
for idx = 1:K
    for noise_level = 2:length(sigma)
        eval("crlb(noise_level) = crlb(noise_level) + CRLB_tdoaOTHR(F, Rb, Ym, Rm, R, beta"+string(idx)+", XYZ, emitter"+string(idx)+", sigma(noise_level), M,index);");
    end
end
err = [];
clf
for idx_SNR = 1:length(sigma)
    filename = ".\real_"+string(sigma(idx_SNR))+".txt";
    [err_1,err_2,err_3,err_4]=textread(filename,'%f%f%f%f','delimiter',',');
    err = [err (mean(err_1) + mean(err_2) + mean(err_3) + mean(err_4))];
end
plot(sigma,crlb,'r^-');
hold on
plot(sigma,1000*err,'ok-','linewidth',1.1,'markerfacecolor',[36, 169, 225]/255);
legend("CRLB",'AMLC');