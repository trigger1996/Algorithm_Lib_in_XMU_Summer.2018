clear;
clc;
clear all;

f = @(x, y) x .* cos(2 * pi * y) + y .* sin(2 * pi * x);	% x属于[-2, 2] y属于[-5, 3]，不对称区间

pc = 0.8;							% 交叉概率
pm = 0.05;							% 变异概率

N = 200;							% 种群上限
L = 21;								% 基因长度
iter = 15;							% 迭代次数，修改迭代次数以查看收敛情况
dna = randi([0, 1], [N, L]);		% 基因

dcd_x = [ 512; 256; 128; 64; 32; 16; 8; 4; 2; 1; 0;    0;   0;   0;   0;  0;  0;  0; 0; 0; 0 ];
dcd_y = [ 0;   0;   0;   0;  0;  0;  0; 0; 0; 0; 1024; 512; 256; 128; 64; 32; 16; 8; 4; 2; 1 ];

for gen = 1 : iter;
	
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
			end

			x1(i, :) = [ dna(i, 1 : d1-1), m(d1 : d2),       dna(i, d2+1 : L) ];
			x2(i, :) = [ m(1 : d1-1),      dna(i, d1 : d2),  m(d2+1 : L) ];
		end
	end
	
	% 变异
	x3 = dna;
	for i = 1: N
		if rand < pm
			x3(i,randi(L)) = randi([0, 1]);
		end
	end
	
	% 计算
	dna = [dna; x1; x2; x3];							% 合并新旧基因
	x = dna * dcd_x * (2 + 2) / 1023 - 2;
	y = dna * dcd_y * (3 + 5) / 2047 - 5;
	fi = f(x, y);										% 计算适应度
	dna = [dna, fi];
	dna = flipud(sortrows(dna, L + 1));					% 对适应度进行排名
	
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

x = dna * dcd_x * (2 + 2) / 1023 - 2;
y = dna * dcd_y * (3 + 5) / 2047 - 5;
ezmesh(f, [-2, 2, -5, 1]); hold on;						% 画出函数图像
plot3(x, y, f(x, y),'ro','linewidth',3)					% 画出最终解的位置
