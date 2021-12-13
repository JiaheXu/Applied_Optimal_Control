function [xs, us, cost,exitflag,output] = trajopt_sqp(xs,us, S, cb, varargin)
% Example code for discrete trajectory optimization using direct collocation with SQP
% @param xs current trajectory guess (n-x-(N+1) matrix)
% @param us current controls guess (c-x-N matrix)
% @param S system properties (dynamics, cost, constraints,etc...)
% @return xs optimized trajectory
% @return us optimized controls
% @return cost computed cost
% @return exitflag -- from fmincon (see fmincon docs)
% @return output -- from fmincon (see fmincon docs)
% Author: Marin Kobilarov, marin(at)jhu.edu

% in previous code
% us = zeros(2,S.N);
% xs = traj(x0, us, S);

    if ~isfield(S, 'f')
        disp('[E] should provide dynamics function f')
        return
    end
    
    if ~isfield(S, 'L')
        disp('[E] should provide cost function L')
        return
    end
    
    % inequality constraint c(t,x,u)<0
    if ~isfield(S, 'con')
        S.con = [];
    end
    
    % final inequality constraint c(x(t_f))<0
    if ~isfield(S, 'conf')
        S.conf = [];
    end
    
    S.n = size(xs,1); % dimension of state
    S.c = size(us,1); % dimension of control
    S.N = size(us,2); % number of control
    
    %size of xs should be n x N+1
    options = optimset('GradObj','on','GradConstr', 'off', ...
    'MaxIter', 10000, 'MaxFunEvals', 20000, 'TolCon',1e-4, 'TolFun', 1e-3, 'TolX', 1e-4);
    
    %params: h, gi, gf, vi, vf, vg, dgi, cb
    x0 = xs(:,1);
    
    % z is the combination of states and controls
    % there should be S.n*S.N parameters for states, anf S.c*S.N parameters for controls
    z = [reshape(xs(:,2:end), S.n*S.N, 1);
         reshape(us, S.c*S.N, 1) ];
    
    [z,cost,exitflag,output] = fmincon(@objfun, z, [], [], [], [], [], [],...
    @nonlcon, options, x0, S,cb, varargin{:})
    
    xs = [x0, reshape(z(1:S.N*S.n), S.n, S.N)];
    us = reshape(z(S.N*S.n+1:end), S.c, S.N);
end

function [f, g] = objfun(z, x0, S, cb, varargin)
% [f,g] is a form of matlab for fmincon objfunc check https://www.mathworks.com/help/optim/ug/fmincon.html
% f is total loss, g is gradient of L (a list)

    xs = [x0, reshape(z(1:S.N*S.n), S.n, S.N)];
    % xs is a S.n * S.N matrix, S.n is the dimension of state

    us = reshape(z(S.N*S.n+1:end), S.c, S.N);
    % us is a S.c * S.N matrix, S.c is the dimension of control

    
    f=0;
    
    for i=1:S.N+1
        if i<S.N+1
            % S.L func is to compute loss in each interval
            [L, Lx, Lxx, Lu, Luu] = S.L(i, xs(:,i), us(:,i),S);
        % control gradients
            uind = S.N*S.n + (i-1)*S.c;
            g(uind+1:uind+S.c) = Lu;
        else
            % S.Lf func is to compute loss in terminal
            [L, Lx, Lxx] = S.Lf(xs(:,i), S);
        end
        f = f + L;
        %set state gradients
        if (i>1)
            %xind is x_index
            xind = (i-2)*S.n;
            g(xind+1:xind+S.n) = Lx;
        end
    end
    
end

function [c, ceq] = nonlcon(z, x0, S, cb, varargin)
% nonlinear constraints
% here we're not exploiting gradient/sparsity info, which is ciritical for efficiency
% c(x) <=0
% ceq(x)=0
% c and ceq are list of value
    xs = [x0, reshape(z(1:S.N*S.n), S.n, S.N)];
    us = reshape(z(S.N*S.n+1:end), S.c, S.N);
    % callback
    if ~isempty(cb)
        cb(xs, us, S);
    end

    ceq = zeros(S.n*S.N, 1);
    c = [];

    for i=1:S.N
        %ind is state_index
        ind = (i-1)*S.n; % discrete dynamics
        ceq(ind+1 : ind+S.n) = xs(:,i+1) - S.f(i, xs(:,i),us(:,i), S);
        if ~isempty(S.con) % state-control inequality constraints
            c = [c; S.con(i, xs(:,i), us(:,i), S)];
        end
    end
% add final inequality constraint
    if ~isempty(S.conf)
        c = [c; S.conf(xs(:,end), S)];
    end
end