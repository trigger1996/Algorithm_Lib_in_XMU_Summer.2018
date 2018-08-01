clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, x取值范围[-1, 2]
%f = @(x) sin(x) + x .* cos(x);		% 函数表达式

pc = 0.9;							% 交叉概率
pm = 0.1;							% 变异概率

N = 125;							% 种群上限
L = 10;								% 基因长度
iter = 150;							% 迭代次数
dna = randi([0, 1], [N, L]);		% 基因，rand()参数1是每个数的取值范围，参数2是矩阵的行列值

dcd = [ 512; 256; 128; 64; 32; 16; 8; 4; 2; 1 ];
%num = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;

% 总迭代
for gen = 1 : iter
	% 交叉
	for i = 1: size(dna, 1)
		if rand < pc							
			d = randi(size(dna, 1));			% 确定另一个交叉的个体
			m = dna(d,:);						% 确定另一个交叉的个体
			d = randi(L-1);						% 确定交叉断点
			x1(i,:) = [dna(i,1:d), m(d+1:L)];	% 新个体 1，这里注意一下，永远是原始种群第i个个体和第随机个个体换        
			x2(i,:) = [m(1:d), dna(i,d+1:L)];	% 新个体 2
		end
	end
	
	% 变异
	x3 = dna;
	for i = 1: size(dna, 1)						% 变异操作
		if rand < pm
			x3(i,randi(L)) = randi([0, 1]);
		end
	end
	
	% 排序和淘汰
	% 轮盘赌法
	dna_t = [dna; x1; x2; x3];							% 合并新旧基因
	fi = f(dna_t * dcd * ((2 - -1) / (2^10 - 1)) - 1);	% 计算适应度
	%dna = [dna, fi];

	% 轮盘赌法
	% 归一化概率值
	norm = 0;
	for i = 1 : size(dna, 1);
		norm = norm + fi(i);
	end
	fi  = fi / norm;
	% 计算累计概率
	pfi(1) = fi(1);
	for i = 2 : size(dna, 1);
		pfi(i) = pfi(i - 1) + fi(i);
	end
	% 转轮盘
	dna = [];
	for i = 1 : N
		t = rand;
		for j = 1 : size(dna_t, 1)
			if t <= pfi(j)								% t < pfi(j)
				dna(i, :) = dna_t(j, :);
				break;
			end
		end
	end
	
end
x = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;			% 对最终种群解码

fplot(f, [-1, 2]); hold on;								% 画出函数图像
plot(x, f(x),'ro','linewidth',3)						% 画出最终解的位置
