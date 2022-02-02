% This code performs the online iteration of the reduced Schwarz method
clear;
% addpath(genpath('../src/'))
addpath('../src/elliptic/')
%% Domain parameters
Lx = 1.0; Nx = 2^(9); dx = Lx/Nx; x0 = 0:dx:Lx; Nx = length(x0);
Ly = 1.0; Ny = 2^(9); dy = Ly/Ny; y0 = 0:dy:Ly; Ny = length(y0);

plotz = 1; 
[xx,yy] = meshgrid(x0,y0);

x_dx = (x0(1:end-1)+x0(2:end))/2; y_dx = y0;
x_dy = x0; y_dy = (y0(1:end-1)+y0(2:end))/2;

%% Domain decomposition
Dx_buffer = 2^(-4); Dx_overlap = 2^(-4); Mx = 4; 
Dy_buffer = 2^(-4); Dy_overlap = 2^(-4); My = 4; 

x_nw_o = 0:Lx/Mx:Lx-Lx/Mx; x_nw_o = max(x_nw_o-Dx_overlap,0);
y_nw_o = 0:Ly/My:Ly-Ly/My; y_nw_o = max(y_nw_o-Dy_overlap,0);
x_se_o = Lx/Mx:Lx/Mx:Lx; x_se_o = min(x_se_o+Dx_overlap,Ly);
y_se_o = Ly/My:Ly/My:Ly; y_se_o = min(y_se_o+Dy_overlap,Ly);


%% Equation parameters
n = 4;
epsilon = 2^(-n);

f = @(u) u.^3;
del_f = @(u) 3*u.^2; 
f_no = 1;

a = @(x,y) 2+sin(2*pi*x).*cos(2*pi*y)...
    +(2+1.8*sin(2*pi*x/epsilon))./(2+1.8*cos(2*pi*y/epsilon))...
    +(2+sin(2*pi*y/epsilon))./(2+1.8*cos(2*pi*x/epsilon));
a_n = 1;

bdy_D_n = @(x) sin(2*pi*x); bdy_no_n = 1;
bdy_D_s = @(x) -sin(2*pi*x); bdy_no_s = 1;
bdy_D_w = @(y) sin(2*pi*y); bdy_no_w = 1;
bdy_D_e = @(y) -sin(2*pi*y);  bdy_no_e = 1;

%%  Dictionary parameters

radius_n = 20;

dim_r = 5;       

N_dict = 64;                                       % number of samples used in dictionary

expmt = 1;



%% Ref soln
dx_ref = 2^(-12);
x_ref = 0:dx_ref:Lx; y_ref = 0:dx_ref:Ly;
Nx_ref = length(x_ref); Ny_ref = length(y_ref);
[xx_ref,yy_ref] = meshgrid(x_ref,y_ref);

load(fullfile('data_elliptic',['u_ref_fno',int2str(f_no),'_eps',int2str(n),...
    '_bdy',int2str(bdy_no_n),int2str(bdy_no_s),int2str(bdy_no_w),int2str(bdy_no_e),...
    '_a',int2str(a_n),'_dx',num2str(dx_ref),'.mat']),'u_ref','t_ref');

figure(110)
mesh(xx_ref,yy_ref,u_ref); zlim([-2,2]);
xlabel('x'); ylabel('y'); zlabel('u');

step = dx/dx_ref;
u_ref_res = u_ref(1:step:end,1:step:end);


%%    Load Dictionary
Dic = cell(Mx,My); Ind = cell(Mx,My); Wet = cell(Mx,My);

for k = 1:My
    for j = 1:Mx
        
        if j~=1 && j~=Mx && k~=1 && k~=My
            
            load(fullfile('data_elliptic',['G_H12_fno',int2str(f_no),...
                '_eps',int2str(n),...
                '_Mx',int2str(Mx),'_My',int2str(My),...
                '_(',int2str(j),',',int2str(k),')',...
                '_Nx',int2str(Nx),'_Ny',int2str(Ny),...
                '_r',num2str(radius_n,2),...
                '_Nsample',int2str(N_dict),...
                '_Nd',int2str(dim_r),...
                '_dxb',num2str(Dx_buffer),'_dxo',num2str(Dx_overlap),...
                '_expmt',int2str(expmt),'.mat']),'u','phi','t_dic');
            
        else
            
            load(fullfile('data_elliptic',['G_H12_fno',int2str(f_no),...
                '_bdynoWESN',int2str(bdy_no_w),int2str(bdy_no_e),int2str(bdy_no_s),int2str(bdy_no_n),...
                '_eps',int2str(n),...
                '_Mx',int2str(Mx),'_My',int2str(My),...
                '_(',int2str(j),',',int2str(k),')',...
                '_Nx',int2str(Nx),'_Ny',int2str(Ny),...
                '_r',num2str(radius_n,2),...
                '_Nsample',int2str(N_dict),...
                '_Nd',int2str(dim_r),...
                '_dxb',num2str(Dx_buffer),'_dxo',num2str(Dx_overlap),...
                '_expmt',int2str(expmt),'.mat']),'u','phi','t_dic');
            
        end
        
        Dic{j,k} = u; Ind{j,k} = phi;

    end
end


knn = 45; % number of neighbors

%% Reduced Schwarz
tic
[u,qq] = reduced_Schwarz_elliptic_l2(x_nw_o,y_nw_o,x_se_o,y_se_o,...
    Dic,Ind,knn,dx);
toc

%%   Error

linf_err_relative = err_inf(u,u_ref_res);
l2_err_relative = err_l2(u,u_ref_res);
h1_err_relative = err_h1(u,u_ref_res,dx);

ax = a(x_dx,y_dx'); ay = a(x_dy,y_dy');
energy_err_relative = err_energy(u,u_ref_res,ax,ay);


%% Plot
figure(11)
mesh(xx,yy,u); 
xlim([0,Lx]); ylim([0,Ly]); zlim([-plotz,plotz]);

