function output = del(tau,del_tau)
    output = [];
    for i = 1:length(del_tau)
        [~,index] = min(abs(tau(i,:) - del_tau(i)));
        output_temp = tau(i,:);
        output_temp(index) = [];
        output(i,:) = output_temp;
    end
end