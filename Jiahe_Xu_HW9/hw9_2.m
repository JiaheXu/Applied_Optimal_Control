clc;
clear;
% Kalman filtering of the double integrator with position measurements
% timing
dt = 1; % time-step
N = 100; % total time-steps
T = N*dt; % final time

% noise terms
sigma_n = 1.5e-5;
sigma_u = 3e-9;
sigma_v = 3e-6;
S.q = diag([sigma_u,sigma_v]); % external disturbance variance

S.r = sigma_n; % measurement noise variance

% F matrix
S.F = [1 -dt;
       0 1];

% G matrix
S.G = [dt;0];

% Q matrix
S.Q = [sigma_v*dt + (1.0/3.0)*sigma_u*(dt^3), -0.5*sigma_u*(dt^2);
    -0.5*sigma_u*(dt^2), sigma_u*dt];

% R matrix
S.R = S.r;

% H
S.H = [1, 0];

% initial estimate of mean and covariance
P = diag([1e-4 1e-12]);
x = [0; 1.7e-7];

xts = zeros(2, N+1); % true states
xs = zeros(2, N+1); % estimated states
Ps = zeros(2, 2, N+1); % estimated covariances
zs = zeros(1, N); % estimated state
pms = zeros(1, N); % measured position

yita_u = sqrt(sigma_u)*randn(1);
betas = yita_u*(1:N); %bias of u

xts(:,1) = x;
xs(:,1) = x;
Ps(:,:,1) = P;

for k=1:N
    yita_v = sqrt(sigma_v)*randn(1);
    u = 0.02 + betas(k) + yita_v;
    xts(:,k+1) = S.F*xts(:,k) + S.G*u; % % true state

    [x,P] = kf_predict(x,P,u,S); % prediction
    z = xts(1,k+1) + sqrt(S.r)*randn; % generate random measurement
    
    [x,P] = kf_correct(x,P,z,S); % correction
    % record result
    xs(:,k+1) = x;
    Ps(:,:,k+1) = P;
    zs(:,k) = z;
end

plot(xts(1,:), '--', 'LineWidth',2)
hold on
plot(xs(1,:), 'g', 'LineWidth',2)
plot(2:N+1,zs(1,:), 'r', 'LineWidth',2)
legend('true', 'estimated','measured')

plot(xs(1,:) + 1.96*reshape(sqrt(Ps(1,1,:)),N+1,1)', '-g')
plot(xs(1,:) - 1.96*reshape(sqrt(Ps(1,1,:)),N+1,1)', '-g')

function [x,P] = kf_predict(x, P, u, S)
    x = S.F*x + S.G*u;
    P = S.F*P*S.F' + S.Q;
end

function [x,P] = kf_correct(x, P, z, S)
    K = P*S.H'*inv(S.H*P*S.H' + S.R);
    P = (eye(length(x)) - K*S.H)*P;
    x = x + K*(z - S.H*x);
end