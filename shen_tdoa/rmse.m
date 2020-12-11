function [rmseIni] = rmse(yE,y,NumMSE,M)
% yE: source estimate; y: true location

AllComb=perms([1:M]); % all possible permutation
rmse=zeros(NumMSE,1);
for idxMSE=1:NumMSE
  for mm=1:M
      rmse(idxMSE)=rmse(idxMSE)+norm(yE(:,AllComb(idxMSE,mm))-y(:,mm))^2;
  end
end
rmseIni=min(rmse);

end