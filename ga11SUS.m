clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, x取值范围[-1, 2]
%f = @(x) sin(x) + x .* cos(x);		% 函数表达式

pc = 0.8;							% 交叉概率
pm = 0.1;							% 变异概率

N = 50;								% 种群上限
L = 20;								% 基因长度
iter = 150;							% 迭代次数
dna = randi([0, 1], [N, L]);		% 基因，rand()参数1是每个数的取值范围，参数2是矩阵的行列值

dcd = [ 524288; 262144; 131072; 65536; 32768; 16384; 8192; 4096; 2048; 1024; 512; 256; 128; 64; 32; 16; 8; 4; 2; 1 ];
%num = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;

% 总迭代
for gen = 1 : iter
	% 交叉
	for i = 1 : N
		if rand < pc
			d  = randi(N);							% 获取1个用于交换的个体
			m  = dna(d, :);

			d1 = randi(L);
			d2 = randi(L);							% 随机得到两个断点
			if d1 > d2
				t = d1;								% 交换变量，这里务必保证断点d1在d2之前
				d1 = d2;
				d2 = t;
				clear t; 
			end
			% 保护，防止出问题
			if d1 >= 20
				d1 = 18;
				d2 = 19;
			end
			if d2 >= 20
				d2 = 19;
			end

			x1(i, :) = [ dna(i, 1 : d1-1), m(d1 : d2),       dna(i, d2+1 : L) ];
			x2(i, :) = [ m(1 : d1-1),      dna(i, d1 : d2),  m(d2+1 : L) ];
		end
	end
	
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
	fi = f(dna_t * dcd * ((2 - -1) / (2^20 - 1)) - 1);	% 计算适应度
	%dna = [dna, fi];

	% 随机均匀采样法SUS
	% 归一化概率值
	fi_min = min(fi);
	fi = fi - fi_min;
	fi = cumsum(fi);									% 累加
	fi = fi / fi(size(fi, 1));							% 累计概率最后一个值一定是最大值，即全部相加，除以它，就归一化了
	
	% 创建指针
	% 指针是均匀隔开的，所以间隔则必然是：
	% 当前种群树 / 目标种群树，即隔多少个取一个
	% 但是对于概率，由于概率是0~1，所以这个值也必须调整，肯定要归一化，可以认为是在1内取多少个，所以
	% 结果是 1 / nSel
	nInd = size(fi, 1);
	nSel = N;
	interval = 1 / nSel;
	susPtr = 0. : interval : 1.;
	susPtr = susPtr + rand / nSel;		% rand / nSel < 1 / nSel，rand / nSel是一个小于区间步长的随机数
	
	% 根据指针取数
	j = 1;
	dna = [];
	for i = 1 : nInd
		if j > nSel
			break;
		end
		if fi(i) > susPtr(j)
			dna(j, :) = dna_t(i, :);
			j = j + 1;
		end
	end
	
end
x = dna * dcd * ((2 - -1) / (2^20 - 1)) - 1;			% 对最终种群解码

fplot(f, [-1, 2]); hold on;								% 画出函数图像
plot(x, f(x),'ro','linewidth',3)						% 画出最终解的位置
