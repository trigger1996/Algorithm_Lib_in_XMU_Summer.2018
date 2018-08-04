clear;
clc;
clear all;

%%
%  上一部分的工作

locationNum = 4;
cityList    = 1 : locationNum;

m = 3;			% 蚂蚁数量
alpha = 1;
beta  = 2;
rho   = 0.5;

% 这里直接用PPT的题目来做，这样这样对不对自己心里比较有数
% https://blog.csdn.net/peiwang245/article/details/78072130
D = [ 0, 3, 1, 2;
	  3, 0, 5, 4; 
	  1, 5, 0, 2;
	  2, 4, 2, 0; ];
  
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
Cnn = Cnn +D(k1, 1);							% 最后要从当前点回到A点，因为复制的矩阵被破坏得差不多了，所以直接用原矩阵算

D0 = D;											% 强迫症

% 这里需要注意，求启发因子的时候需要求倒数，所以这边必须处理掉所有的0
D(find(D == 0)) = 1e-6;
Eta = 1 ./ D;									% 启发因子，取距离的导数，这个是用来计算转移概率的的

% 初始化信息素矩阵
tau = ones(locationNum) * (m / Cnn);			% 信息素矩阵对角线为0，其余值为m/Cnn
tau = tau - eye(locationNum) * (m / Cnn);		% 这里用全1矩阵减去对角矩阵得到

%%
% 这一部分的工作，放置蚂蚁，计算概率

antTraj = zeros(m, locationNum);			% 蚂蚁轨迹
antTraj(:, 1) = randi([1, 4], [1, 3]);		% 设置蚂蚁初始位置，蚂蚁开始时是可以放在同一个城市的

for i = 1 :locationNum - 1											% 轮次，本来原计划外环算蚂蚁的，然后发现问题太多了，还是先算轮次
	for j = 1 : m
		allowList = cityList(~ismember(cityList, antTraj(j, :)));	% 这行是查找哪些城市能去的
																	% ismemeber是查询antTraj中哪些成员属于cityList，返回下标
																	% 然后取反以后就是哪些不是，不在轨迹上的就是可以去的
		Pt = [];	
		for k = 1 : size(allowList, 2)
			currentLoc = antTraj(j, i);
			nextLoc    = allowList(k);
			Pt(k, :) = tau(currentLoc, nextLoc)^alpha + Eta(currentLoc, nextLoc)^beta;
																	% 根据公式计算转移概率，Pt的顺序和allowList是一样的
		end
		Pt = Pt ./ sum(Pt);											% 根据公式计算概率
		% 轮盘赌求下一个移动的城市
		Pr = cumsum(Pt);											% 求累计概率
		seq = find(rand <= Pr);
		seq = min(seq);												% 计算轮盘指针，轮盘指针必须刚好比累计概率小，累计概率略大一点的就是指针位置
		antTraj(j, i + 1)  = allowList(seq);						% 轨迹写入
	end
end
antTraj = [ antTraj, antTraj(:, 1) ];								% 轨迹加上初始点，蚂蚁必须回到初始点