clc;
clear; 
N = 100;
x0 = [10; 0];
% x0 = [10; 5];
% x0 = [10; -5];
Pf = eye(2);
R = 0.01;
dt = 0.1;
% x_dot = Ax + Bu + w
A = [ 1        , dt ; 
      0.2 * dt , 1 - 0.5*dt ];

B = [0 ; 1];
w = [0 ; 0.1];

% Value function params
P = cell(N);
P{N} = Pf;

b = cell(N);
b{N} = [0;0];

c = cell(N);
c{N} = 0;

% trajectory 
x = zeros(2,N);
x(:,1) = x0;

% control law
K = cell(N);
k = cell(N);
u = zeros(1,N);

% get P_i, b_i, c_i, K_i, k_i
for i = N-1:-1:1
    % every variable has the same notatoin within Q2
    % A B w R are constants here
    S = ( R + B' * P{i + 1} * B);
    inv_S = inv( R + B' * P{i + 1} * B);

    K_i = - inv_S * B' * P{i + 1} * A;
    k_i = - inv_S * B' * (P{i + 1} * w + b{i + 1});
    % determine P_i, b_i, c_i
    % there is no Q_i in Q3
    P_i = A'*( P{i+1} - P{i+1}' * B * inv_S * B' * P{i+1}  )* A;
    b_i = (( w' * P{i+1} + b{i+1}' )*( A - B*inv_S*(B')*P{i+1}*A ))';
    M = (B'*( P{i+1}*w + b{i+1}));

    c_i = 1/2 * w'*P{i+1}*w - 1/2*(M')*inv_S*M +(b{i+1}')*w + c{i+1};

    K{i} = K_i;
    k{i} = k_i;
    
    P{i} = P_i;
    b{i} = b_i;
    c{i} = c_i;
    
end

for i = 1:N-1
    % get x_i to integrate to x_i+1
    x_i = x(:,i);  % i-th column
    
    % get the control law
    K_i = K{i};
    k_i = k{i};
    u_i = K_i * x_i + k_i;
    % get x_i+1
    x_next = A * x_i + B * u_i + w;
    
    x(:,i+1) = x_next;
    u(i) = u_i;
end

% Plotting results
f_traj = figure(1);
plot(x(1,:), x(2,:), 'DisplayName', 'trajectory'); 
title(sprintf('Trajectory: R = %.3f', R));
xlabel('p'); ylabel('v');
legend();
hold on;
plot(x0(1), x0(2), 'r*', 'DisplayName', 'start'); 
grid on;

f_ctrl = figure(2);
plot(u);
title(sprintf('Control: R = %.3f', R));
xlabel('i'); ylabel('u');
grid on;




