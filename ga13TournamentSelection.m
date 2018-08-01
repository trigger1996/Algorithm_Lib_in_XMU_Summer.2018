clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, x取值范围[-1, 2]
%f = @(x) sin(x) + x .* cos(x);		% 函数表达式

pc = 0.8;							% 交叉概率
pm = 0.05;							% 变异概率

N = 50;								% 种群上限
L = 10;								% 基因长度
iter = 6;							% 迭代次数
dna = randi([0, 1], [N, L]);		% 基因，rand()参数1是每个数的取值范围，参数2是矩阵的行列值

dcd = [ 512; 256; 128; 64; 32; 16; 8; 4; 2; 1 ];
%num = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;

% 总迭代
for gen = 1 : iter
	% 交叉
	for i = 1: N
		if rand < pc							
			d = randi(N);						% 确定另一个交叉的个体
			m = dna(d,:);						% 确定另一个交叉的个体
			d = randi(L-1);						% 确定交叉断点
			x1(i,:) = [dna(i,1:d), m(d+1:L)];	% 新个体 1，这里注意一下，永远是原始种群第i个个体和第随机个个体换        
			x2(i,:) = [m(1:d), dna(i,d+1:L)];	% 新个体 2
		end
	end
	
	% 变异
	x3 = dna;
	for i = 1: N								% 变异操作
		if rand < pm
			x3(i,randi(L)) = randi([0, 1]);
		end
	end
	
	% 排序和淘汰
	% 竞标赛法选择
	% 竞标赛法就是一次随机选出几个，进行比赛，然后适应度最高的塞进去
	% 而且允许重复选择
	% 写起来很容易
	dna_t = [dna; x1; x2; x3];							% 合并新旧基因
	fi = f(dna_t * dcd * ((2 - -1) / (2^10 - 1)) - 1);	% 计算适应度
	dna_t = [dna_t, fi];
	
	nIden = size(dna_t, 1);
	nSel  = N;
	dna   = [];
	for i = 1 : nSel
		tournSel = randi([1, nIden], [1, 5]);				% 每次随机挑选出来多少个，进行竞赛
		tournSel = sort(tournSel);							% 排序适应度就是比赛了
		
		compti = dna_t(tournSel, :);						% 选中的序号挑出
		compti = flipud(sortrows(compti, L + 1));			% 对适应度进行排名
		dna(i, :) = compti(1, :);							% 最好的一个被选出
		
	end
	dna = dna(:, 1:L);	
end
x = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;			% 对最终种群解码

fplot(f, [-1, 2]); hold on;								% 画出函数图像
plot(x, f(x),'ro','linewidth',3)						% 画出最终解的位置
