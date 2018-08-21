clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, x取值范围[-1, 2]
x_max = 2;
x_min  = -1;

% 初始化粒子参数
iter = 200;

m = 25;
v_max = 0.5;
v_min = -v_max;

c1 = 1;
c2 = 2;

% 初始化粒子群
x(:, 1) = randi([x_min * 1000, x_max * 1000], [1, m]) ./ 1000;
v(:, 1) = randi([v_min * 1000, v_max * 1000], [1, m]) ./ 1000;

% 个体最优值
pbest(:, 1) = x(:, 1);              % 这里注意一下，pbest比较的是因变量f(x)，但是记录的是自变量x的值
% 种群最优值
gbest = x(1, 1);
for i = 1 : m
    if f(gbest) < f(pbest(i, 1))
        gbest = pbest(i, 1);
    end
end

for gen = 1 : iter

    % 更新每个粒子发现的最好位置
    for i = 1 : m
        if f(x(i, 1)) > f(pbest(i, 1))
            pbest(i, 1) = x(i, 1);
        end
    end
    
    % 更新整个群体目前发现的最好的位置
    for i = 1 : m
        if f(gbest) < f(pbest(i, 1))
            gbest = pbest(i, 1);
        end
    end
    
    % 更新速度
    v(:, 1) = v(:, 1) + c1 * rand * (pbest(:, 1) - x(:, 1)) + c2 * rand * (gbest - x(:, 1));
    
    % 速度限制
   v(find(v(:, 1) > v_max)) = v_max;
   v(find(v(:, 1) < v_min)) = v_min;   
    
    
    % 更新位置
    x(:, 1) = x(:, 1) + v(:, 1);
    
    % 位置限制
   x(find(x(:, 1) > x_max)) = x_max;
   x(find(x(:, 1) < x_min)) = x_min;       
    
end


figure(1);
fplot(f, [-1, 2]); hold on;								% 画出函数图像
plot(x, f(x), 'ro', 'linewidth', 3);                    % 作出粒子的最终位置
plot(gbest, f(gbest), 'bo', 'linewidth', 3);            % 作出最好的粒子的历史位置
