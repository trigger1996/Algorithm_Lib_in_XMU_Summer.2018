clear;
clc;
clear all;

%
% 在地图构建完毕的基础上
% 进行地图编码
%
locationNum = 10;
pos = randi([0, 10000], [locationNum, 2]) / 100.;			% 产生一个0~100的(x, y)坐标向量，randi的参数前者是范围，是矩阵的Size

figure;
for i = 1 : locationNum
	plot(pos(i, 1), pos(i, 2),'ro','linewidth',3); hold on;	% 画出最终解的位置
	text(pos(i, 1) + 1, pos(i, 2) + 1, num2str(i), 'FontSize', 12);
end

costmap = zeros(10, 10);
for i = 1 : locationNum
	for j = 1 : locationNum
		dx = pos(i, 1) - pos(j, 1);
		dy = pos(i, 2) - pos(j, 2);
		costmap(i, j) = sqrt(dx*dx + dy*dy);
	end
end

%% 遗传算法部分
% 产生初始随机数列数列
N = 50;
L = locationNum;

pc = 0.8;							% 交叉概率
pm = 0.25;							% 变异概率

iter = 200;					% 迭代次数
for i = 1 : N
	dna(i, :) = randperm(10);								% 产生不重复的随机数列
end

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
	for i = 1: N								% 变异操作，这边的变异不能rand了，因为序号不能重复，所以宁可两个位置相互交换
		if rand < pm
			k1 = randi([1, L]);
			k2 = randi([1, L]);
			t = x3(i, k1);
			x3(i, k1) = x3(1, k2);
			x3(i, k2) = t;
		end
	end
	
	
	% 排序和淘汰
	dna_t = [dna; x1; x2; x3];							% 合并新旧基因
	
	% 计算适应度
	% 其实计算没那么复杂，就是查表
	dna_t = [dna, dna(:, 1) ];							% 最后必须回到开始的地方
	for i = 1 : N
		fi(i, 1) = 0;
		for j = 2 : L + 1
			k1 = dna_t(i, j);
			k2 = dna_t(i, j - 1);
			fi(i, 1) = fi(i, 1) + costmap(k1, k2);
		end
	end
	
	dna = [dna, fi];
	dna = sortrows(dna, L + 1);							% 对适应度进行排名，注意这回排序是从低到高
	
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

route = [ dna(1, :), dna(1, 1) ];
for i = 2 : L + 1
	posK1 = route(:, i - 1);
	posK2 = route(:, i);
	
	posX = pos(posK1, 1);
	posY = pos(posK1, 2);
	posU = pos(posK2, 1);
	posV = pos(posK2, 2);
	
	quiver(posX, posY, posU - posX, posV - posY);
	
end

