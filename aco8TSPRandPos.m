clear;
clc;
clear all;

%%
%  上一部分的工作

iter_max = 100;

locationNum = 8;
cityList    = 1 : locationNum;

m = ceil(locationNum * 1.5);								% 蚂蚁数量，取城市数量的1.5倍整数
alpha = 1;
beta  = 2;
rho   = 0.3;

% 这次这个改成随机坐标
pos = randi([0, 10000], [locationNum, 2]) / 1000.;			% 产生一个0~1000的(x, y)坐标向量，randi的参数前者是范围，是矩阵的Size
% 计算邻接矩阵
D = zeros(locationNum);
for i = 1 : locationNum
	for j = 1 : locationNum
		dx = pos(i, 1) - pos(j, 1);
		dy = pos(i, 2) - pos(j, 2);
		D(i, j) = sqrt(dx*dx + dy*dy);
	end
end

  
% 计算初始信息素
% 贪心算法算最短路径作为参考
k1 = 1;				% 开始从A点出发
k2 = 1;
Cnn = 0;

D0 = D;
D0(find(D == 0)) = inf;
% 运行贪心算法，其实这一部分人工手算也是可以的，毕竟不是重点
for i = 1 : size(D0, 1) - 1
	Cnn = Cnn + min(D0(k1, :));
	k2 = find(D0(k1, :) == min(D0(k1, :)));		% 找到这行最小值的元素的下标——目标是找到贪心算法算出来的下一个位置，下一个位置要从这里开始计算
	D0(:, k1) = inf;							% 上回走过的点不能去了，所以整列设成无限大，直接禁止
	k1 = k2;									% 传递参数，准备下一轮计算
end
Cnn = Cnn + D(k1, 1);							% 最后要从当前点回到A点，因为复制的矩阵被破坏得差不多了，所以直接用原矩阵算

D0 = D;											% 强迫症

% 这里需要注意，求启发因子的时候需要求倒数，所以这边必须处理掉所有的0
D(find(D == 0)) = 1e10;
Eta = 1 ./ D;									% 启发因子，取距离的导数，这个是用来计算转移概率的的

% 初始化信息素矩阵
tau = ones(locationNum) * (m / Cnn);			% 信息素矩阵对角线为0，其余值为m/Cnn
tau = tau - eye(locationNum) * (m / Cnn);		% 这里用全1矩阵减去对角矩阵得到

%%

for iter = 1 : iter_max

	antTraj = zeros(m, locationNum);									% 蚂蚁轨迹
	antTraj(:, 1) = randi([1, locationNum], [1, m]);					% 设置蚂蚁初始位置，蚂蚁开始时是可以放在同一个城市的

	for i = 1 :locationNum - 1											% 轮次，本来原计划外环算蚂蚁的，然后发现问题太多了，还是先算轮次
		for j = 1 : m
			allowList = cityList(~ismember(cityList, antTraj(j, :)));	% 这行是查找哪些城市能去的

			Pt = [];													% 加个保护
			for k = 1 : size(allowList, 2)
				currentLoc = antTraj(j, i);
				nextLoc    = allowList(k);
				Pt(k, :) = tau(currentLoc, nextLoc)^alpha + Eta(currentLoc, nextLoc)^beta;
																		% 根据公式计算转移概率，Pt的顺序和allowList是一样的
			end
			Pt = Pt ./ sum(Pt);											% 根据公式计算概率
			% 采用排名法选择下一次的路径
			% 原来这里是轮盘赌，然后一直不收敛，所以试一下排名
			Pr = [ allowList', Pt ];
			Pr = flipud(sortrows(Pr, 2));
			for k = 1: size(Pr, 1)
				if rand < Pr(k, 2)
					antTraj(j, i + 1)  = Pr(k, 1);
					break;
				end
			end
			if antTraj(j, i + 1) == 0
				antTraj(j, i + 1) = Pr(1, 1);
			end
		end
	end
	antTraj = [ antTraj, antTraj(:, 1) ];								% 轨迹加上初始点，蚂蚁必须回到初始点



	% 根据邻接矩阵计算总路程
	Cxx = zeros(m, 1);
	for i = 1 : m
		for j = 1 : locationNum
			currentLoc = antTraj(i, j);
			nextLoc    = antTraj(i, j + 1);
			Cxx(i, :) = Cxx(i, :) + D(currentLoc, nextLoc);				% 很多代码上面两个变量是不写的，所以看起来贼难看懂，但是其实都是根据邻接矩阵算距离
		end
	end

	% 更新信息素矩阵
	% 首先算蒸发，蒸发好算
	tau = (1 - rho) * tau;

	% 然后根据蚂蚁的路程来算信息素矩阵
	% 根据公式，每次更新都是直接加就好了本次蚂蚁运动的总路程的倒数就好了
	for i = 1 : m
		for j = 1 : locationNum
			currentLoc = antTraj(i, j);
			nextLoc    = antTraj(i, j + 1);
			tau(currentLoc, nextLoc) = tau(currentLoc, nextLoc) + 1 / Cxx(i, :);
		end
	end

	% 记录最佳路径
	seq_out = find(min(Cxx));
	bestTrajCxx(iter, :) = Cxx(seq_out, :);
	
end

% 最后输出那只路径最短的蚂蚁的路径即可
seq_out = find(min(Cxx));
bestTraj = antTraj(seq_out, :);

figure(1);
for i = 1 : locationNum
	plot(pos(i, 1), pos(i, 2),'ro','linewidth',3); hold on;						% 画出坐标点
	text(pos(i, 1) + .2, pos(i, 2) + .2, num2str(i), 'FontSize', 12);				% 标出坐标序号
end

% 作出路线
for i = 2 : locationNum + 1
	posK1 = bestTraj(:, i - 1);
	posK2 = bestTraj(:, i);
	
	posX = pos(posK1, 1);
	posY = pos(posK1, 2);
	posU = pos(posK2, 1);
	posV = pos(posK2, 2);
	
	quiver(posX, posY, posU - posX, posV - posY);
	
end

% 作出收敛曲线
figure(2);
plot(1 : iter, 	bestTrajCxx);

% 输出最佳路径
bestTraj
