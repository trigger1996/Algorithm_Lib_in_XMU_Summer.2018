clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, xȡֵ��Χ[-1, 2]
x_max = 2;
x_min  = -1;

% ��ʼ�����Ӳ���
iter = 200;

m = 25;
v_max = 0.5;
v_min = -v_max;

c1 = 1;
c2 = 2;

% ��ʼ������Ⱥ
x(:, 1) = randi([x_min * 1000, x_max * 1000], [1, m]) ./ 1000;
v(:, 1) = randi([v_min * 1000, v_max * 1000], [1, m]) ./ 1000;

% ��������ֵ
pbest(:, 1) = x(:, 1);              % ����ע��һ�£�pbest�Ƚϵ��������f(x)�����Ǽ�¼�����Ա���x��ֵ
% ��Ⱥ����ֵ
gbest = x(1, 1);
for i = 1 : m
    if f(gbest) < f(pbest(i, 1))
        gbest = pbest(i, 1);
    end
end

for gen = 1 : iter

    % ����ÿ�����ӷ��ֵ����λ��
    for i = 1 : m
        if f(x(i, 1)) > f(pbest(i, 1))
            pbest(i, 1) = x(i, 1);
        end
    end
    
    % ��������Ⱥ��Ŀǰ���ֵ���õ�λ��
    for i = 1 : m
        if f(gbest) < f(pbest(i, 1))
            gbest = pbest(i, 1);
        end
    end
    
    % �����ٶ�
    v(:, 1) = v(:, 1) + c1 * rand * (pbest(:, 1) - x(:, 1)) + c2 * rand * (gbest - x(:, 1));
    
    % �ٶ�����
   v(find(v(:, 1) > v_max)) = v_max;
   v(find(v(:, 1) < v_min)) = v_min;   
    
    
    % ����λ��
    x(:, 1) = x(:, 1) + v(:, 1);
    
    % λ������
   x(find(x(:, 1) > x_max)) = x_max;
   x(find(x(:, 1) < x_min)) = x_min;       
    
end


figure(1);
fplot(f, [-1, 2]); hold on;								% ��������ͼ��
plot(x, f(x), 'ro', 'linewidth', 3);                    % �������ӵ�����λ��
plot(gbest, f(gbest), 'bo', 'linewidth', 3);            % ������õ����ӵ���ʷλ��
