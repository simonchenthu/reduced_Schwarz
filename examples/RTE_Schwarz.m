% This code uses classical Schwarz iteration to solve RTE
clear;
% addpath(genpath('../src/'))
addpath('../src/RTE/')
%% Domain decomposition
L = 3.0; Nx = 3*2^(11); dx = L/Nx; x0 = 0:dx:L; Nx = length(x0);
t = [0,0.25,0.75,1.25,1.75,2.25,2.75]; s = [0.25,0.75,1.25,1.75,2.25,2.75,3]; 
dx_overlap = 2^(-3);

t = [t(1),t(2:end)-dx_overlap];          % left end points of patches
s = [s(1:end-1)+dx_overlap,s(end)];      % right end points of patches

Mx = length(t);

%% Equation and solver
n = 4;
epsilon=2^(-n);

Nv = 2^7; [v0,w0] = legendre_quad(Nv-1); v0 = v0'; w0 = w0'; 
w_bdy = [w0';1;1]; 
[xx,vv] = meshgrid(x0,v0); 
ploty = 4; plotz = 40;

% sigma = @(x) 2+sin(2*pi*x); sigma_n = 1;
sigma = @(x) ones(size(x)); sigma_n = 0;
sigma_x = sigma(x0)';

% v<0, x = L
bdy_vn = 2+sin(2*pi*v0(1:Nv/2)'); bdy_xL = 3; bdy_l_no = 1;    
% bdy_vn = 2-(v0(1:Nv/2).^2)'; bdy_xL = 3; bdy_l_no = 2;
% bdy_vn = 3-(v0(1:Nv/2).^2)'; bdy_xL = 1; bdy_l_no = 3;

% v>0, x=0
bdy_vp = 3+sin(2*pi*v0(Nv/2+1:end)'); bdy_x0 = 2; bdy_r_no = 1;  
% bdy_vp = 3-(v0(Nv/2+1:end).^2)'; bdy_x0 = 1; bdy_r_no = 2;
% bdy_vp = 1.5-(v0(Nv/2+1:end).^2)'; bdy_x0 = 2.5; bdy_r_no = 3;
                                                                      

%% Ref soln
dx_ref = 2^(-14); Nv_ref = 2^10; x_ref = 0:dx_ref:L;
[v_ref,w_v] = legendre_quad(Nv_ref-1); w_v = w_v'; 
w_x = [dx_ref/2,dx_ref*ones(1,L/dx_ref-1),dx_ref/2]';
[xx_ref,vv_ref] = meshgrid(x_ref,v_ref);
load(fullfile('data_RTE',['G_ref_eps',int2str(n),'_bdy',int2str(bdy_l_no),int2str(bdy_r_no),...
    '_sigma',int2str(sigma_n),'_Lx',num2str(L,2),'_dx',num2str(dx_ref),'_Nv',int2str(Nv_ref),...
    '.mat']),'f_ref','theta_ref');

theta_bdy = [bdy_x0;bdy_xL]; %[left;right]
f_bdy = [bdy_vn;bdy_vp]; %[v<0;v>0]

%% Schwarz iteration
tic

[f,theta,res,q] = Schwarz_RTE(epsilon,sigma_x,f_bdy,theta_bdy,t,s,dx,v0,w0);

t_s = toc


%%   Error  %%%
f_intp = interp2(xx,vv,f,xx_ref,vv_ref,'spline');
theta_intp = interp1(x0,theta,x_ref,'spline');
l2_error_relative = sqrt(w_v*(f_ref-f_intp).^2*w_x + (theta_ref-theta_intp).^2*w_x)...
    /sqrt(w_v*f_ref.^2*w_x + theta_ref.^2*w_x);

%% Plot
figure(3)
subplot(2,1,1); 
plot(x0,theta); 
ylim([0,ploty]);
title('Schwarz Solution');

subplot(2,1,2);
mesh(xx,vv,f); 
zlim([0,plotz]);

%% Save
save(fullfile('data_RTE',['G_Schwarz_eps',int2str(n),...
    '_bdy',int2str(bdy_l_no),int2str(bdy_r_no),...
    '_sigma',int2str(sigma_n),'_Lx',num2str(L,2),'_dx',num2str(dx),'_Nv',int2str(Nv),...
    '_dx_overlap',num2str(dx_overlap),'.mat']),'f','theta','t_s','res','q','l2_error_relative');

