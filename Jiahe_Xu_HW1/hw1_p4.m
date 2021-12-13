%%%%%%%%%%%%
% hw1_p4
%%%%%%%%%%%%
clc;
clear;

syms x1 x2 real;
x = [x1 x2]';
% a
f = (1-x1)^2 + 200*(x2 - x1*x1 )^2;
x0 = [0,0]' ;
% figure;
% ezsurfc(f,[-4 4]);
% hold off;
% 1 stands for gradient decent;
method = 1;
% last 1 stands for constant stepsize
% last 2 stands for variable stepsize
%a
disp(['############################################### a']);
optimize(x0,x,f,method , 1 ); %figure1
optimize(x0,x,f,method , 2 ); %figure2


%b
disp(['############################################### b']);
method = 2;
optimize(x0,x,f,method , 1 ); %figure3
optimize(x0,x,f,method , 2 ); %figure4


% c
disp(['############################################### c']);
f = x1 * exp( -x1*x1 - 0.5*x2*x2 ) + x1*x1 / 10 + x2*x2 / 10;
x0 = [1.5 , 3]';
% figure;
% ezsurfc(f,[-4 4]);
% hold off;
%2 stands for newton method
method = 1;
optimize(x0, x, f ,method,2); %figure5
method = 2;
optimize(x0, x, f ,method,2); %figure6

function optimize(x0, x, f, method , var_stepsize )

    
    figure
    %ezsurfc(f,[-2 2],100)
    %hold on

    grad = jacobian(f,x).';
    hess = hessian(f,x);
    
    fx = matlabFunction(f,'vars',{x});
    gx = matlabFunction(grad,'vars',{x});
    hx = matlabFunction(hess,'vars',{x});
    
    x=x0;
    step = .0001;
    xs = x;
    Steps = 0;
        
    if(method == 1)
        while norm(gx(x)) > 1e-8
            x = gradient_decent(fx, gx, x, step , var_stepsize );
            xs = [xs, x];
            Steps = Steps+1;
        end
        plot(xs(1,:), xs(2,:), '-ob');
        if(var_stepsize == 1)
            disp(['gradient decent with fixed stepsize:']);
        else
            disp(['gradient decent with variable stepsize:']);
        end
        disp(['x=[' num2str(x') '] f=' num2str(fx(x)) ]);
        disp(['Steps=' num2str(Steps)]);
        hold off;
    else
        while norm(gx(x)) > 1e-8  && Steps < 1000000
            x = newton_method(fx, gx, hx, x , step , var_stepsize);
            xs = [xs, x];
            Steps = Steps+1;
        end
        plot(xs(1,:), xs(2,:), '-ob');
         if(var_stepsize == 1)
            disp(['newton method with fixed stepsize:']);
        else
            disp(['newton method with variable stepsize:']);
        end

        disp(['x=[' num2str(x') '] f=' num2str(fx(x)) ]);
        disp(['Steps=' num2str(Steps)]);
        hold off;
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of optimize()

function x = gradient_decent(fx, gx, x, stepsize, var_stepsize)

    g = gx(x);
    if norm(g) < 1e-8
        return;
    end
    d = -g;
    if(var_stepsize == 1)
        a = stepsize;
    else
        a =  armijo(fx, x, g, d);
    end
    x = x + a*d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of gradient_decent()

function x = newton_method(fx, gx, hx, x , stepsize, var_stepsize)

    g = gx(x);
    h = hx(x);

    if norm(g) < 1e-8
        return;
    end

    d = -inv(h)*g;
    if(var_stepsize == 1)
        a = stepsize;
    else
        a =  armijo(fx, x, g, d);
    end
    x = x + a*d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% newton_method

function a = armijo(fx, x, g, d)
    sigma = .01;
    beta = .25;
    s = 1;
    f = fx(x);
    for m=0:1000
        a = beta^m*s;
        if (f - fx(x + a*d) > -sigma*(beta^m)*s*g'*d)
            return;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of armijo