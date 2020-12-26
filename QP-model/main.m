clear all
R = 6371.2;
noise = 500;
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
M = 4;
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
XYZ = zeros(M,3);
%Bei Jing
[x0 y0 z0] = LGLTtoXYZ(116.41,39.90,R);
XYZ(1,:) = [x0 y0 z0];
%Wu Han
[x0 y0 z0] = LGLTtoXYZ(114.31,30.59,R);
XYZ(2,:) = [x0 y0 z0];
%Shang Hai
[x0 y0 z0] = LGLTtoXYZ(121.47,31.23,R);
XYZ(3,:) = [x0 y0 z0];
%Tokyo
[x0 y0 z0] = LGLTtoXYZ(139.69,35.69,R);
XYZ(4,:) = [x0 y0 z0];
% %Seoul
% [x0 y0 z0] = LGLTtoXYZ(126.58,37.33,R);
% XYZ(5,:) = [x0 y0 z0];
% %Qing Dao
% [x0 y0 z0] = LGLTtoXYZ(120.3826,36.0671,R);
% XYZ(6,:) = [x0 y0 z0];                  

%Jia Yi
[x0 y0 z0] = LGLTtoXYZ(120.4491,23.4801,R);
emitter1 = [x0 y0 z0]';
beta = zeros(1,M);
x = emitter1';
for k = 1:M
    beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
end
beta1 = beta;
emitter(:,1) = emitter1;
% beta2 = [0.130743022187950,0.391197329272201,0.520620592610200,0.0738255619888921,0.168515720771200];
%Yin Chuan
[x0 y0 z0] = LGLTtoXYZ(106.2309,38.4872,R);
emitter2 = [x0 y0 z0]';
beta = zeros(1,M);
x = emitter2';
for k = 1:M
    beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
end
beta2 = beta;
emitter(:,2) = emitter2;
% beta1 = [0.486164591596301,0.310478405869350,0.177057112583560,0.00974700715358300,0.143452972295470];
%Qiqi Haer
[x0 y0 z0] = LGLTtoXYZ(123.9182,47.3543,R);
emitter3 = [x0 y0 z0]';
% beta3 = [0.373296850690061,0.105761408290660,0.140327847255840,0.134834308670000,0.315574458947491];
beta = zeros(1,M);
x = emitter3';
for k = 1:M
    beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
end
beta3 = beta;
emitter(:,3) = emitter3;
%Hong Kong
[x0 y0 z0] = LGLTtoXYZ(114.16,22.28,R);
emitter4 = [x0 y0 z0]';
beta = zeros(1,M);
x = emitter4';
for k = 1:M
    beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
end
beta4 = beta;
emitter(:,4) = emitter4;
K = 4;

tau = [];
for i = 1:K
    sigma_t = noise * 10^-9 * 3 * 10^5 ;
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
    eval("param.P"+string(i)+" = permutation(K);");
end
tau_0 = tau(1:M-1,:);
P_tau = [];
for i = 1:M-1
    eval("P_tau = [P_tau;tau(i,:)*param.P"+string(i)+"];");
end
%Initialize the linear cut sets
ini = '';
for i = 1:M-1
     ini = ini + "param.cut"+string(i)+"(:,:,1) = zeros(K,K);";
end
eval(ini);
cvx_solver GUROBI_2
x_rec = [];
for iter = 1:100
    if K ~= 1
        [P_tau0,param] = IP(M,G,t,P_tau,K,param);
    end
    obj = trace((G(1:M-1,1:M)*t - P_tau0)'*inv_Omega*(G(1:M-1,1:M)*t - P_tau0));
    fprintf("iter:%d obj:%2.8f K:%d\n",iter,obj,K);
    [t_sum,obj_sum,location] = solve_GPGD(M,N,F,Rb,Rm,Ym,P_F,R,P_Rb,P_Rm,P_Ym,G,P_tau0,inv_Omega,upper,max_dis,min_dis,XYZ,plt,K);    
	t = t_sum;
    index = find(obj_sum >= 1e1);
    dif_index = find(obj_sum <= 1e0);
    if isempty(index)
        %Record the location of the last iterartion and break
        for i = 1:length(dif_index)
            eval("x_rec = [x_rec;location.x"+string(dif_index(i))+"];");
        end
        break
    end
    if ~isempty(dif_index)
        %Update P_tau and K
        for i = 1:length(dif_index)
            P_tau = del(P_tau,P_tau0(:,dif_index(i)));
        end
        K = size(P_tau,2);
        %Record the location
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
%Compute the localization error
[P] = compute_err(x_rec,emitter);
x = P*x_rec;
for i = 1:size(emitter,2)
    err(i) = norm(x(i,:) - emitter(:,i)');
end
%Save the results
fid=fopen("real_M4_"+string(noise)+".txt","a+");
fprintf(fid,"%2.4f",err(1));
for i = 2:size(emitter,2)
    fprintf(fid,",%2.4f",err(i));
end
fprintf(fid,"\n");
fclose(fid);
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

