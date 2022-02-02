%%% Finite difference solver with Anderson acceleration for outer layer and
%%% monotone iteration for semielliptic eqn
function [f,theta] = nonlinearRTE_FD_aa_monotone(f_bdy,theta_bdy,epsilon,sigma_x,...
                                                                x,w0,v0,dx)
Nx = length(x);
Nv = length(v0);
v0 = v0';

% coefficients in exponential finite difference
temp11 = exp(dx./(epsilon*v0(1:Nv/2))); 
temp12 = 1-exp(dx./(epsilon*v0(1:Nv/2)));

temp21 = exp(-dx./(epsilon*v0(Nv/2+1:end))); 
temp22 = 1-exp(-dx./(epsilon*v0(Nv/2+1:end)));

% constant used in monotone iteration
gamma = max(nthroot(f_bdy(1:Nv/2)./sigma_x(end),4));
gamma = max(gamma,max(nthroot(f_bdy(Nv/2+1:end)./sigma_x(1),4)));
gamma = max(gamma,max(theta_bdy)); 
lambda = 6*gamma^3;

% elliptic operator in monotone iteration
elliptic_op = spdiags([-ones(Nx-2,1), (2+lambda*dx^2/epsilon^2)*ones(Nx-2,1), -ones(Nx-2,1)],-1:1,Nx-2,Nx-2);




% initial f
f = ones(Nv,Nx); f_hat = w0/2*f;
theta = compute_theta(theta_bdy,f_hat,epsilon,sigma_x,dx,Nx,elliptic_op,lambda);

TOL_i = 2e-3; err_i = 1;         % smaller TOL_i may stablize the acceleration
while err_i > TOL_i
    
    theta_new = fixed_pt_map(theta,sigma_x,f_bdy,theta_bdy,epsilon,dx,Nx,Nv,w0,elliptic_op,lambda...
                                              ,temp11,temp12,temp21,temp22);
    err_i = norm(theta-theta_new,2)/norm(theta,2);
    theta = theta_new;
end

it_p = 10;                        % suitable m may stablize the acceleration; depends on the dimension of the solution manifold
beta = 1; r = 0.001;

theta_acc = zeros(Nx,it_p+1);
theta_acc(:,1) = theta;
for i = 2:it_p+2
    
    theta = fixed_pt_map(theta,sigma_x,f_bdy,theta_bdy,epsilon,dx,Nx,Nv,w0,elliptic_op,lambda...
                                              ,temp11,temp12,temp21,temp22);
    theta_acc(:,i) = theta;
    
end

Gc = theta_acc(2:end-1,2:end) - theta_acc(2:end-1,1:end-1);
Fc = Gc(:,2:end) - Gc(:,1:end-1); Gc = Gc(:,2:end);

f1 = Gc(:,end);
g1 = theta_acc(2:end-1,end);

% err = norm(theta_acc(:,end)-theta_acc(:,end-1),2)/norm(theta_acc(:,end),2);

TOL = 1e-8; i = 1; err = 1;

while err > TOL
    
    i_update = i - it_p*floor((i-1)/it_p);
    
    g2 = fixed_pt_map(theta,sigma_x,f_bdy,theta_bdy,epsilon,dx,Nx,Nv,w0,elliptic_op,lambda...
                                              ,temp11,temp12,temp21,temp22);
    g2 = g2(2:end-1);
    f2 = g2 - theta(2:end-1);
    
    Fc(:,i_update) = f2-f1; Gc(:,i_update) = g2-g1;
    
    gamma = Fc\f2;
    theta_new = g2 - Gc*gamma;                         % undamped
    
    theta_new = theta_new - (1-beta)*(f2-Fc*gamma);    % damping
    
    err_aa = norm(theta(2:end-1)-theta_new,2)/norm(theta(2:end-1),2);
    
    f1 = f2; g1 = g2;
    if err < err_aa*r || err_aa > 0.1
        err = norm(theta(2:end-1)-g2,2)/norm(theta(2:end-1),2);
        theta(2:end-1) = g2;
    else
        err = err_aa;
        theta(2:end-1) = theta_new; 
    end
    i = i+1; 
    
end

theta4 = sigma_x.*theta.^4;
theta4_mid = (theta4(1:end-1)+theta4(2:end))/2;

f = transport(theta4_mid,f_bdy,Nx,Nv,temp11,temp12,temp21,temp22);
theta = theta';

end

function theta = fixed_pt_map(theta,sigma_x,f_bdy,theta_bdy,epsilon,dx,Nx,Nv,w0,elliptic_op,lambda...
                                              ,temp11,temp12,temp21,temp22)
% UpdatingF
theta4 = sigma_x.*theta.^4;
theta4_mid = (theta4(1:end-1)+theta4(2:end))/2;
f = transport(theta4_mid,f_bdy,Nx,Nv,temp11,temp12,temp21,temp22);

% UpdatingTheta;
f_hat = w0/2*f;
theta = compute_theta(theta_bdy,f_hat,epsilon,sigma_x,dx,Nx,elliptic_op,lambda);
    
end

function theta = compute_theta(theta_bdy,f_hat,epsilon,sigma_x,dx,Nx,elliptic_op,lambda)

% Calculate theta given rho by monotone iteration

%lower solution
theta = zeros(Nx-2,1);
%initialization
TOL = 1e-8; err = 1;
while err > TOL
    
    theta_new = monotonic_map(theta,theta_bdy,f_hat,sigma_x,elliptic_op,dx,epsilon,lambda);
    
    err = norm(theta-theta_new,2)/norm(theta_new,2);
    
    theta = theta_new;
    
end

theta = [theta_bdy(1);theta;theta_bdy(2)];

end

function theta = monotonic_map(theta,theta_bdy,f_hat,sigma_x,elliptic_op,dx,epsilon,lambda)

d = lambda*theta - sigma_x(2:end-1).*theta.^4 + f_hat(2:end-1)'; d = dx^2/epsilon^2*d;
d(1) = d(1) + theta_bdy(1); d(end) = d(end) + theta_bdy(2);
theta = elliptic_op\d;

end

function f = transport(theta4,f_bdy,Nx,Nv,temp11,temp12,temp21,temp22)

f=zeros(Nv,Nx);

f(1:Nv/2,end)=f_bdy(1:Nv/2);
for i=Nx-1:-1:1
    f(1:Nv/2,i) = temp11.*f(1:Nv/2,i+1) + temp12.*theta4(i);
end

f(Nv/2+1:end,1)=f_bdy(Nv/2+1:end);
for i=2:Nx
    f(Nv/2+1:end,i) = temp21.*f(Nv/2+1:end,i-1) + temp22.*theta4(i-1);
end

end