clc;
clear;

syms tf v
a = sqrt(6)/2;

F = @(x) [ x(2) * ( 4*( (exp(a*x(1))-exp(-a*x(1)))/(exp(a*x(1))+exp(-a*x(1))) )^2 -(tan(a*x(1)))^2 -5) + 10 ;
    x(2)/(2*a)*(-8* (exp(a*x(1))-exp(-a*x(1)))/(exp(a*x(1))+exp(-a*x(1)))-2*tan(a*x(1)) )+5*x(1) - 15];
all_ans = [];
num_ans = 0;

for i=0:0.1:10
    for j=-10:0.1:10
        x0 = [i;j];
        x = fsolve(F,x0);
        x= x';
        res = F(x);
        if( res'*res < 1e-3)
            found = 0;
            for k=1:num_ans
                if( ( x - all_ans(k,:))' * (x-all_ans(k,:)) < 1e-3 )
                    found = 1;
                end
            end
            if(found == 0)
                all_ans = [all_ans;x];
                num_ans = num_ans+1;
            end
        end
    end
end
size(all_ans)
all_ans