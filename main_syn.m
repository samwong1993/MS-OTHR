clear all
R = 6371.2;
plt = 0;
if plt == 1
    figure('color','k')
    ha=axes('units','normalized','position',[0 0 1 1]);
    uistack(ha,'down')
    img=imread('background.jpg');
    image(img)
    colormap gray
    set(ha,'handlevisibility','off','visible','off');
    rotate3d on
    hold on
    %Add legend
    %emitter
    point1 = scatter3(0,0,0,50,'filled','r');
    %sensors
    point2 = scatter3(0,0,0,'filled','c');
    %initial point
    point3 = scatter3(0,0,0,50,'filled','b');
    %recover location
    point4 = scatter3(0,0,0,50,'k*');
    %generated sequence
    point5 = scatter3(0,0,0,5,'m');
    %plot earth
    npanels = 72;
    alpha   = 1;
    image_file = 'earth.jpg';
    [x0, y0, z0] = ellipsoid(0, 0, 0, R, R, R);
    globe = surf(x0, y0, -z0, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    cdata = imread(image_file);
    set(gca, 'NextPlot','add', 'Visible','off');
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    axis off;
    axis equal;
    axis auto;
end
M = 6;
K = 4;
N = M*(M-1)/2;
Omega = covariance(1,M);
inv_Omega = inv(Omega(1:M-1,1:M-1));
delta = 0;
Rm = 6650;
Ym = 100;
Rb = Rm - Ym;
fc = 10;
f = 15;
F = f/fc;
%perturbation on ionosphere
P_Rm = Rm;
P_Ym = Ym;
P_Rb = P_Rm - P_Ym;
P_fc = fc;
P_f = f;
P_F = P_f/P_fc;
[max_dis,min_dis,upper] = beta_bound(M,F,R,Rb,Rm,Ym);
[G] = generate_G(N,M);
[emitter,XYZ,beta0] = generator(M,K,F,R,Rb,Rm,Ym,max_dis,min_dis);

for k = 1:K
    eval("emitter"+string(k)+" = emitter(:,"+string(k)+");");
    eval("beta"+string(k)+" = beta0("+string(k)+",:);");
end
      

tau = [];
for i = 1:K
    sigma_t = 1000 * 10^-9 * 3 * 10^5 ;
    noise_t0 = randn(M,1);
    noise_t = (sigma_t*G*noise_t0)';
    eval("tau"+string(i)+" = generate_tau(M,F,R,Rb,Rm,Ym,emitter"+string(i)+",XYZ) + noise_t;");
    tau = [tau,eval("tau"+string(i)+"'")];
end

t_0 = [];
for i = 1:K
    [A B C] = ABC(F,R,Rb,Rm,Ym,eval("beta"+string(i)));
    eval("[t"+string(i)+" D] = PD(A,B,C,beta"+string(i)+",R,Rb);");
    t_0 = [t_0,eval("t"+string(i)+"'")];
end

t = zeros(M,K);
for i = 1:M-1
    eval("P_"+string(i)+" = permutation(K);");
end

tau_0 = tau(1:M-1,:);
P_tau = [tau(1,:)*P_1;tau(2,:)*P_2;tau(3,:)*P_3;tau(4,:)*P_4;tau(5,:)*P_5];

ini = '';
for i = 1:M-1
     ini = ini + "param.cut"+string(i)+"(:,:,1) = zeros(K,K);";
end
eval(ini);
cvx_solver mosek
x_rec = [];
for iter = 1:100
    [P_tau0 param] = IP(M,G,t,P_tau,K,param);
    obj = trace((G(1:M-1,1:M)*t - P_tau0)'*inv_Omega*(G(1:M-1,1:M)*t - P_tau0));
    fprintf("obj:%2.8f K:%d\n",obj,K);
    diag1 = trace((G(1:M-1,1:M)*t - P_tau0)'*inv_Omega*(G(1:M-1,1:M)*t - P_tau0));
    [t_sum,obj_sum,location] = solve_GPGD(M,N,F,Rb,Rm,Ym,P_F,R,P_Rb,P_Rm,P_Ym,G,P_tau0,inv_Omega,upper,max_dis,min_dis,XYZ,plt,K);
    diag2 = trace((G(1:M-1,1:M)*t_sum - P_tau0)'*inv_Omega*(G(1:M-1,1:M)*t_sum - P_tau0));
    
	t = t_sum;
    index = find(obj_sum >= 1e1);
    dif_index = find(obj_sum <= 1e1);
    if isempty(index)
        for i = 1:length(dif_index)
            eval("x_rec = [x_rec;location.x"+string(dif_index(i))+"];");
        end
        break
    end
    P_tau = P_tau0(:,index);
	K_old = K;
    K = size(P_tau,2);
    if K ~= K_old
        for i = 1:length(dif_index)
            eval("x_rec = [x_rec;location.x"+string(dif_index(i))+"];");
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

[P] = compute_err(x_rec,emitter);
x = P*x_rec;
for i = 1:size(emitter,2)
    err(i) = norm(x(i,:) - emitter(:,i)');
end
fprintf("(%2.2f",err(1))
for i = 2:size(emitter,2)
    fprintf(",%2.2f",err(i))
end
fprintf(")\n")


if plt == 1
    scatter3(emitter1(1),emitter1(2),emitter1(3),50,'filled','r')
    scatter3(emitter2(1),emitter2(2),emitter2(3),50,'filled','r')
    scatter3(emitter3(1),emitter3(2),emitter3(3),50,'filled','r')
    scatter3(emitter4(1),emitter4(2),emitter4(3),50,'filled','r')
    scatter3(x(1),x(2),x(3),50,'k*')
    %text(emitter(1),emitter(2),emitter(3),'e')
end
for i = 1:M
    if plt == 1
        scatter3(XYZ(i,1),XYZ(i,2),XYZ(i,3),'filled','c')
        temp = strcat('s ',num2str(i));
        text(XYZ(i,1),XYZ(i,2),XYZ(i,3),temp);
    end
    fprintf("Sensor %d:(%2.2f,%2.2f,%2.2f)\n",i,XYZ(i,1),XYZ(i,2),XYZ(i,3))
end

