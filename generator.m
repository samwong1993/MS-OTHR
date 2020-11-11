%% Generate emitter and sensors
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [emitter,XYZ,beta0] = generator(M,K,F,R,Rb,Rm,Ym,max_dis,min_dis)
emitter(1,:) = randn(1,3);
emitter(1,:) = R*emitter(1,:) / norm(emitter(1,:));
for i = 1:M
    while(1)
        XYZ(i,:) = randn(1,3);
        XYZ(i,:) = R*XYZ(i,:) / norm(XYZ(i,:));
        if norm(emitter(1,:) - XYZ(i,:))<max_dis&norm(emitter(1,:) - XYZ(i,:))>min_dis
            break
        end
    end
end




for k = 2:K
    while(1)
        index = zeros(1,M);
        emitter(k,:) = randn(1,3);
        emitter(k,:) = R*emitter(k,:) / norm(emitter(k,:));
        for i = 1:M
            if norm(emitter(k,:) - XYZ(i,:))<max_dis&norm(emitter(k,:) - XYZ(i,:))>min_dis
                index(i) = 1;
            else
                break
            end
        end
        if all(index) == 1
            break
        end
    end
end

%Calculate corresponing beta of x
for k = 1:K
while(1)
beta = zeros(1,M);
x = emitter(k,:);
for i = 1:M
    beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,i);
end
[A B C] = ABC(F,R,Rb,Rm,Ym,beta);
[P D] = PD(A,B,C,beta,R,Rb);
for i = 1:M
    penalty(i) = norm(XYZ(i,:)-x,2)^2 - 4*R^2*(1 - cos(D(i)/R))/2;
end
if sum(penalty) < 1e-5
    beta0(k,:) = beta;
    break;
end
end
end
emitter = emitter';
end