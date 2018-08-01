clear;
clc;
clear all;

f = @(x) sin(10 * pi * x) ./ x;
%f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, x取值范围[0.5, 3]
%f = @(x) sin(x) + x .* cos(x);		% 函数表达式

pc = 0.8;							% 交叉概率
pm = 0.05;							% 变异概率

N = 50;								% 种群上限
L = 12;								% 基因长度
iter = 15;							% 迭代次数，修改迭代次数以查看收敛情况
dna = randi([0, 1], [N, L]);		% 基因

dcd = [ 2048; 1024; 512; 256; 128; 64; 32; 16; 8; 4; 2; 1 ];
%num = dna * dcd * ((2 - -1) / (2^L - 1)) - 1;

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
	dna = [dna; x1; x2; x3];							% 合并新旧基因
	fi = f(dna * dcd * ((3 - 0.5) / (2^L - 1)) + 0.5);	% 计算适应度
	dna = [dna, fi];
	dna = flipud(sortrows(dna, L + 1));					% 对适应度进行排名
	%dna = sortrows(dna, L + 1);						% 修改1：如果改变排序，则变成最小值
	
	% 自然选择
	% 如果子代个数太多了，就开始淘汰
	while size(dna, 1) > N
		d = randi(size(dna, 1));
		if rand < (d - 1) / size(dna, 1)				% 排名法，越靠后的淘汰几率越高
			dna(d,:) = [];
			fi(d, :) = [];
		end
	end
	dna = dna(:, 1:L);
end
x = dna * dcd * ((3 - 0.5) / (2^L - 1)) + 0.5;			% 对最终种群解码

fplot(f, [0.5, 3]); hold on;							% 画出函数图像
plot(x, f(x),'ro','linewidth',3)						% 画出最终解的位置
