clc;
clear;

l1 = 2;
l2 = 2;
po = 1/2*[2;2];
r = 1/2;

steps = 200;
theta1s = linspace(0, 2*pi, steps);
theta2s = linspace(0, 2*pi, steps);
valid_tips = [];

for i = 1:steps
    for j = 1:steps
        th1 = theta1s(i);
        th2 = theta2s(j);
        tip = [cos(th1)*l1+cos(th1+th2)*l2-po(1);
               sin(th1)*l1+sin(th1+th2)*l2-po(2)];
        if (r - norm(tip) <= 0 )
            valid_tips = [ valid_tips , [th1;th2]];
        end
    end
end
%valid_tips= valid_tips';
figure;
title('problem3 (a)');
valid_tips = rad2deg(valid_tips);
plot(valid_tips(1,:), valid_tips(2,:), '.', 'MarkerSize', 14);
xlabel('theta1 (deg)'); ylabel('theta2 (deg)');
xlim([0, 360]); ylim([0, 360]);
valid_tips = [];
valid_map = zeros(steps,steps);

for i = 1:steps
    for j = 1:steps
        th1 = theta1s(i);
        th2 = theta2s(j);

        d1 = [cos(th1);sin(th1)];
        s1 = min(max(0,po'*d1),l1);
        c1 = s1*d1;

        A2 = l1*d1;
        B2 = [cos(th1)*l1+cos(th1+th2)*l2;
               sin(th1)*l1+sin(th1+th2)*l2];
        d2 = [cos(th1+th2);sin(th1+th2)];
        s2 = (po-A2)'*d2;
        s2 = min(max(0,s2),l2);
        c2 = A2+s2*d2;
        if ( norm(c1-po) >= r && norm(c2-po) >= r)
            valid_tips = [ valid_tips , [th1;th2]];
        end
    end
end
%valid_tips= valid_tips';
figure;
title('problem3 (b)');
valid_tips = rad2deg(valid_tips);
plot(valid_tips(1,:), valid_tips(2,:), '.', 'MarkerSize', 14);
xlabel('theta1 (deg)'); ylabel('theta2 (deg)');
xlim([0, 360]); ylim([0, 360]);