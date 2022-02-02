% Online iteration of the reduced order Schwarz method
clear;
% addpath(genpath('../src/'))
addpath('../src/RTE/')
%% Domain parameters
L = 3.0; Nx = 3*2^(11); dx = L/Nx; x0 = 0:dx:L; Nx = length(x0);
Nv = 2^7; [v0,w0] = legendre_quad(Nv-1); v0 = v0'; w0 = w0'; w_bdy = [w0';1;1]; 
[xx,vv] = meshgrid(x0,v0); 
ploty = 4; plotz = 40;

%% Domain decomposition
t = [0,0.25,0.75,1.25,1.75,2.25,2.75]; s = [0.25,0.75,1.25,1.75,2.25,2.75,3]; 
dx_overlap = 2^(-3); dx_buffer = 2^(-2);

t = [t(1),t(2:end)-dx_overlap];          % left end points of patches
s = [s(1:end-1)+dx_overlap,s(end)];      % right end points of patches

Mx = length(t);


%% Equation parameters
n = 4;
epsilon=2^(-n);  % Knudsen number

% sigma = @(x) 2+sin(2*pi*x); sigma_n = 1;
sigma = @(x) ones(size(x)); sigma_n = 0;  % medium

% boundary condition

%%% v<0, x = L %%%
bdy_vn = 2+sin(2*pi*v0(1:Nv/2)'); bdy_xL = 3; bdy_r_no = 1;    
% bdy_vn = 2-(v0(1:Nv/2).^2)'; bdy_xL = 3; bdy_r_no = 2;
% bdy_vn = 3-(v0(1:Nv/2).^2)'; bdy_xL = 1; bdy_r_no = 3;

%%% v>0, x=0 %%%
bdy_vp = 3+sin(2*pi*v0(Nv/2+1:end)'); bdy_x0 = 2; bdy_l_no = 1;  
% bdy_vp = 3-(v0(Nv/2+1:end).^2)'; bdy_x0 = 1; bdy_l_no = 2;
% bdy_vp = 1.5-(v0(Nv/2+1:end).^2)'; bdy_x0 = 2.5; bdy_l_no = 3;

%% Dictionary parameters
N_sample = 64;
radius = 25;
expmt = 2;

%% Ref soln
% dx_ref = 2^(-14); Nv_ref = 2^10; x_ref = 0:dx_ref:L;
% [v_ref,w_v] = legendre_quad(Nv_ref-1); w_v = w_v'; 
% w_x = [dx_ref/2,dx_ref*ones(1,L/dx_ref-1),dx_ref/2]';
% [xx_ref,vv_ref] = meshgrid(x_ref,v_ref);
% load(fullfile('data_RTE',['G_ref_eps',int2str(n),'_bdy',int2str(bdy_l_no),int2str(bdy_r_no),...
%     '_sigma',int2str(sigma_n),'_Lx',num2str(L,2),'_dx',num2str(dx_ref),'_Nv',int2str(Nv_ref),...
%     '.mat']),'f_ref','theta_ref');
% 
% figure(11)
% subplot(2,1,1); 
% plot(x_ref,theta_ref); 
% ylim([0,ploty]); 
% title('Reference Solution'); 
% xlabel('$x$','Interpreter','latex'); 
% ylabel('$T$','Interpreter','latex');
% 
% subplot(2,1,2); 
% mesh(xx_ref,vv_ref,f_ref); 
% zlim([0,plotz]); 
% xlabel('$x$','Interpreter','latex'); 
% ylabel('$v$','Interpreter','latex'); 
% zlabel('$I$','Interpreter','latex');


%%    Load Dictionary
Dic = cell(2,Mx); Ind = zeros(Nv+2,N_sample,Mx);

% left patch
load(fullfile('data_RTE',['G_bdy',int2str(bdy_l_no),...
    '_eps',int2str(n),'_sigma',int2str(sigma_n),...
    '_Lx',num2str(L,2),'_dx',num2str(dx),'_Nv',int2str(Nv),...
    '_Mx',int2str(Mx),'(',int2str(1),')',...
    '_r',num2str(radius,2),'_Nsample',int2str(N_sample),...
    '_dx_overlap',num2str(dx_overlap),'_dx_buffer',num2str(dx_buffer),...
    '_expmt',int2str(expmt),'.mat']),'f_dic','theta_dic','phi_dic');
Dic{1,1} = f_dic; Dic{2,1} = theta_dic; Ind(:,:,1) = phi_dic;

% right patch
load(fullfile('data_RTE',['G_bdy',int2str(bdy_r_no),...
    '_eps',int2str(n),'_sigma',int2str(sigma_n),...
    '_Lx',num2str(L,2),'_dx',num2str(dx),'_Nv',int2str(Nv),...
    '_Mx',int2str(Mx),'(',int2str(Mx),')',...
    '_r',num2str(radius,2),'_Nsample',int2str(N_sample),...
    '_dx_overlap',num2str(dx_overlap),'_dx_buffer',num2str(dx_buffer),...
    '_expmt',int2str(expmt),'.mat']),'f_dic','theta_dic','phi_dic');
Dic{1,Mx} = f_dic; Dic{2,Mx} = theta_dic; Ind(:,:,Mx) = phi_dic;

% interior patches
for j = 2:Mx-1

    load(fullfile('data_RTE',['G_eps',int2str(n),'_sigma',int2str(sigma_n),...
        '_Lx',num2str(L,2),'_dx',num2str(dx),'_Nv',int2str(Nv),...
        '_Mx',int2str(Mx),'(',int2str(j),')',...
        '_r',num2str(radius,2),'_Nsample',int2str(N_sample),...
        '_dx_overlap',num2str(dx_overlap),'_dx_buffer',num2str(dx_buffer),...
        '_expmt',int2str(expmt),'.mat']),'f_dic','theta_dic','phi_dic');
    
    Dic{1,j} = f_dic; Dic{2,j} = theta_dic; Ind(:,:,j) = phi_dic;

end


knn = 5;

tic
%% Reduced Schwarz
[f,theta,qq,discrep,iter_res] = reduced_Schwarz_RTE(t,s,Dic,Ind,knn,dx,w_bdy);
toc

%%   Error
% f_intp = interp2(xx,vv,f,xx_ref,vv_ref,'spline');
% theta_intp = interp1(x0,theta,x_ref,'spline');
% l2_error_relative = sqrt(w_v*(f_ref-f_intp).^2*w_x + (theta_ref-theta_intp).^2*w_x)...
%     /sqrt(w_v*f_ref.^2*w_x + theta_ref.^2*w_x);

%% Plot    
figure(14)

subplot(2,1,1); plot(x0,theta); ylim([0,ploty]);
title('Patched Approximate Solution'); 
xlabel('$x$','Interpreter','latex'); 
ylabel('$T$','Interpreter','latex');

subplot(2,1,2); 
mesh(xx,vv,f); 
zlim([0,plotz]); 
xlabel('$x$','Interpreter','latex'); 
ylabel('$v$','Interpreter','latex'); 
zlabel('$I$','Interpreter','latex');



