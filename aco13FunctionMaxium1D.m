clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, x取值范围[-1, 2]

interval_max = 2;
interval_min = -1;

iter_max = 1000;		% 迭代次数

m = 50;					% 蚂蚁数量
alpha = 5;				% 信息素的幂因子，这个是用来算转移概率的
%beta  = 1;				% 路程的幂因子，这个是算转移概率的，但是连续函数不打算用这个东西
rho = 0.9;				% 信息素挥发因子

%dx  = (interval_max - interval_min) * 0.1;
						% 局部搜索区间，这个变量后面被弃用了

%
% 这里作两个变动

% 转移概率中，路程作为常数，不使用
% 同时，信息素越强烈的地方，转移概率越低，tau用本轮最大值减去当前值，再归一化
% 所以P = (tau_max - tau) / sum(tau);
% 后面发现这么算P会很小，不太对
% 所以P = (tau_max - tau) / tau_max;

% 信息素的计算，启发因子直接使用函数值本身
% tau = (1 - rho) * tau + F(x);
% 蚂蚁的移动时随机的，如果如果转移概率小，则在自身附近移动，否则在全局随机移动

% 设置蚂蚁的初始位置
x = rand(1, m) * (interval_max - interval_min) + interval_min;		% 在 [区间最小值, 区间最大值] 范围上生成随机数，因为rand是0~1，所以可以这么做
x_init = x;															% 记录初始值，画出初始解

% 初始化信息素矩阵
tau = f(x);															% 初始化信息素就是f(x)

for iter = 1 : iter_max
	% 计算转移概率
	tau_max = max(tau);
	P = (tau_max - tau).^alpha / tau_max;							% (tau_max - tau) / sum(tau)，后来发现这样取值太小了，所以干脆除以最大值
	
	% 蚂蚁开始跑
	for i = 1 : m
		if rand < P(i)												% 满足转移概率则全局搜索
			pos_neg = sign(rand-0.5);								% 网络上的办法，产生随机的正负号
			x(i) = x(i) + pos_neg * (interval_max - interval_min) * rand;
		else														% 否则局部搜索
			pos_neg = sign(rand-0.5);								% 网络上的办法，产生随机的正负号
			x(i) = x(i) + pos_neg * 1 / iter * rand;				% x = x + pos_neg * dx * rand，这个地方会很大程度影响算法收敛性
																	% x = x + pos_neg * 1 / iter * rand
																	% 原文采用后者，收敛性很好，如果采用前者，就完全不收敛
		end
	end
	
	for i = 1 : m
		if x(i) < interval_min
			x(i) = interval_min;
		end
		if x(i) > interval_max
			x(i) = interval_max;
		end
	end
	
	% 更新信息素
	tau = (1 - rho) * tau +  f(x);
	
end

figure(1); 
fplot(f, [-1, 2]); hold on
plot(x_init, f(x_init),'ko','linewidth',2)								% 画出初始解的位置

figure(2);
fplot(f, [-1, 2]); hold on
plot(x, f(x),'ro','linewidth',3)									% 画出最终解的位置
