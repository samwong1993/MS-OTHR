function output = del(tau,del_tau)

output = [];
for i = 1:length(del_tau)
    [~,index] = min(abs(tau(i,:) - del_tau(i)));
    output(i,:) = setdiff(tau(i,:),tau(i,index));
end

end
