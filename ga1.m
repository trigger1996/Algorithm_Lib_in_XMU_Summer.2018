    clear
    clc
    close all
    
	f = @(x) x .* sin(10 * pi * x)
    %f = @(x) sin(x) + x .* cos(x);  % 函数表达式
    ezplot(f, [0, 2*pi])             % 画出函数图像
     
    N = 50;                          % 种群上限
    ger = 100;                       % 迭代次数
    L = 5;                           % 基因长度
    pc = 0.8;                        % 交叉概率
    pm = 0.1;                        % 变异概率
    dco = [10000; 1000; 100; 10 ;1]; % 解码器
    dna = randi([0, 9], [N, L]);     % 基因
    hold on
    x = dna * dco / 99999 * 2 * pi;  % 对初始种群解码
    plot(x, f(x),'ko','linewidth',3) % 画出初始解的位置
     
    x1 = zeros(N, L);                % 初始化子代基因，提速用
    x2 = x1;                         % 同上
    x3 = x1;                         % 同上
    fi = zeros(N, 1);                % 初始化适应度，提速
     
    for epoch = 1: ger               % 进化代数为100
        for i = 1: N                 % 交叉操作
            if rand < pc
               d = randi(N);            % 确定另一个交叉的个体
               m = dna(d,:);            % 确定另一个交叉的个体
               d = randi(L-1);          % 确定交叉断点
               x1(i,:) = [dna(i,1:d), m(d+1:L)];  % 新个体 1        
               x2(i,:) = [m(1:d), dna(i,d+1:L)];  % 新个体 2
            end
        end
        x3 = dna;
        for i = 1: N                           % 变异操作
            if rand < pm
                x3(i,randi(L)) = randi([0, 9]);
            end
        end
        dna = [dna; x1; x2; x3];               % 合并新旧基因
        fi = f(dna * dco / 99999 * 2 * pi);    % 计算适应度，容易理解
        dna = [dna, fi];
        dna = flipud(sortrows(dna, L + 1));    % 对适应度进行排名
        while size(dna, 1) > N                 % 自然选择
            d = randi(size(dna, 1));           % 排名法
            if rand < (d - 1) / size(dna, 1)
                dna(d,:) = [];
                fi(d, :) = [];
            end
        end
        dna = dna(:, 1:L);
    end
    x = dna * dco / 99999 * 2 * pi;            % 对最终种群解码
    plot(x, f(x),'ro','linewidth',3)           % 画出最终解的位置
    disp(['最优解为x=',num2str(x(1))]);
    disp(['最优值为y=',num2str(fi(1))]);