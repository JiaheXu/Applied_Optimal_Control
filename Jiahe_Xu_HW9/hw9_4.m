clc;
clear;

s = 10;
% true shape parameter (i.e. a symmetric cup)
x_true = [1; 1; 0; 0; 0; 0];

% prior parameters:
%m = [1.2; 1.3; 1; 1; 1; 1];
%m = [1.2; 1.2; 1.2; 1.2; 1.2; 1.2];
%m = [1.3; 1.3; 1.3; 1.3; 1.3; 1.3];
m = [1.2; 1.3; 1; 1; 1; 1]*2
%m = [1.2; 1.3; 1; 1; 1; 1]*10;

P0 = diag([16,16,16,16,16,16]);

% plot true
gt = ezsurf(@(p,q)shape(p, q, x_true),[-s,s]);
alpha(gt, 0.3)
hold on
% measurement standard dev
std = 20;

% #of measurements
k = 8;
% generate random measurements
p = 4*s*(rand(k,1) - .5);
q = 4*s*(rand(k,1) - .5);
z = shape(p, q, x_true) + randn(k,1)*std;
% estimate optimal parameters x
R = diag(repmat(std^2, k, 1));

H = shape_basis(p, q);
x = m;
P = P0;
for i = 1:k/2
    start = 2*i-1;
    goal = 2*i;

    R_i = R( start : goal , start : goal );
    H_i = H(start : goal,:);
    Z_i = z(start : goal);
    P = inv(inv(P) + (H_i')*inv(R_i)*H_i);
    Ki = (P*H_i')*inv(R_i);
    x = x + Ki*(Z_i - H_i*x);
end
x
% plot estimated

ge = ezsurf(@(p,q)shape(p,q,x),[-s,s]);
alpha(ge, .8)

function f = shape_basis(p, q)
% quadratic function, although could be any shape
    f = [p.^2, q.^2, p.*q, p, q, ones(size(p))];
end
    
function z = shape(p, q, x)
    z = shape_basis(p, q)*x;
end