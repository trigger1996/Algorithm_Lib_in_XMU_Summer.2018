clear;
clc;
clear all;

%
% 在地图构建完毕的基础上
% 进行地图编码
%
locationNum = 10;
pos = randi([0, 10000], [locationNum, 2]) / 100.;			% 产生一个0~100的(x, y)坐标向量，randi的参数前者是范围，是矩阵的Size

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
N = 25;
L = locationNum;

pc = 0.6;							% 交叉概率
pm = 0.05;							% 变异概率

iter = 200;					% 迭代次数
for i = 1 : N
	dna(i, :) = randperm(10);								% 产生不重复的随机数列
end

for gen = 1 : iter
	
	% 交叉
	for i = 1: N
		d1 = randi(N);
		d2 = randi(N);
		% 这里为了偷懒，直接植入了别的函数的swap，详见intercross.m
		a = dna(d1, :);
		b = dna(d2, :);
		
		if ~isempty(find(a(1, :) == 0))
			a = randperm(10);						% 有时候会出现全0的数组，这里的任务就是，如果找到全0的数组，则重新赋值
		end
		if ~isempty(find(b(1, :) == 0))				% 本来可以直接赋值dna变量里的，但是这边这么做，就是主要怕掉到局部最优里
			b = randperm(10);
		end		
		
		if rand < pc							
			r1=randsrc(1,1,[1:L]);
			r2=randsrc(1,1,[1:L]);
			if r1~=r2
				a0=a;b0=b;							% 总的交换之前保存一个
				s=min([r1,r2]);
				e=max([r1,r2]);
				for j=s:e
					a1=a;b1=b;						% 每一轮交换之前保存一个
					a(j)=b0(j);						% 单个交换
					b(j)=a0(j);
					x=find(a==a(j));				% 单个交换完看一下，当前序列中，有没有和交换完的重复的
					y=find(b==b(j));
					j1=x(x~=j);						% 必然会find找到自己，但是除了自己，如果还有别的相同的，则用i1记录下来，
					j2=y(y~=j);						% 而且每次只交换一个，所以至多2个相同的，一个是自己，一个是原有的，所以如果只有自己，那么i1为空，同理，i2也是这样
					if ~isempty(j1)					% 如果i1非空
						a(j1)=a1(j);				% 一模一样的数字有两个，一个是刚替换的，一个是原来就有的，那么就取修改那个原来就有的，原来就有的改成前者，即刚替换的之前的数字，这样形成交换
					end
					if ~isempty(j2)
						b(j2)=b1(j);
					end
				end
			end
		end
		x1(i, :) = a;
		x2(i, :) = b;
	end
	
	% 变异
	x3 = dna;
	for i = 1: N								% 变异操作，这边的变异不能rand了，因为序号不能重复，所以宁可两个位置相互交换
		if rand < pm
			k1 = randi([1, L]);
			k2 = randi([1, L]);
			t = x3(i, k1);
			x3(i, k1) = x3(i, k2);
			x3(i, k2) = t;
		end
	end
	
	
	% 排序和淘汰
	dna = [dna; x1; x2; x3];							% 合并新旧基因
	
	% 计算适应度
	% 其实计算没那么复杂，就是查表
	dna = [dna, dna(:, 1) ];							% 最后必须回到开始的地方
	for i = 1 : N * 4
		fi(i, 1) = 0;
		for j = 2 : L + 1
			k1 = dna(i, j);
			k2 = dna(i, j - 1);
			fi(i, 1) = fi(i, 1) + costmap(k1, k2);
		end
		fi(i, 1) = 1 / fi(i, 1);						% 因为这边是求最小，所以这边必须取倒数
	end
	
	%dna = [dna, fi];
	
	% SUS法选择
	%fi_min = min(fi);
	pfi = fi;											% 归一化概率，这边减去一个最小值是为了避免受到负数或者过大的正数的影响, pfi = fi - fi_min
	pfi = cumsum(pfi);									% 后面发现因为上面求了倒数，所以如果在收敛的时候，如果都减去一个最小值
	pfi = pfi / pfi(size(pfi, 1));
	dna_t = [ dna, pfi ];								% 趁着还没排序，绑定编码和适应度函数，这样在后面好操作
	
	nIdn = size(fi, 1);									% 产生指针
	nSel = N;
	interval = 1 / nSel;
	susPtr = 0. : interval : 1.;
	susPtr = susPtr + rand / nSel;						% rand / nSel < 1 / nSel，rand / nSel是一个小于区间步长的随机数
		
	if ~isempty(find(isnan(pfi)))
		gen = gen + 1;
		gen = gen - 1;
	end
	
	% 根据指针取数
	j = 1;
	dna = [];
	for i = 1 : nIdn
		if j > nSel
			break;
		end
		if dna_t(i, L + 2) > susPtr(j)
			dna(j, :) = dna_t(i, 1 : L);
			j = j + 1;
		end
	end
	
	% 记录fi的最优值以查看收敛情况
	fi_best(gen) = fi(1, :);
		
end

%% 作图

% 作出位置
figure;
for i = 1 : locationNum
	plot(pos(i, 1), pos(i, 2),'ro','linewidth',3); hold on;	% 画出最终解的位置
	text(pos(i, 1) + 1, pos(i, 2) + 1, num2str(i), 'FontSize', 12);
end

% 作出路线
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

% 作出收敛情况
figure;
plot(1 : 200, fi_best);

