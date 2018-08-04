clear;
clc;
clear all;

locationNum = 4;
m = 3;			% 蚂蚁数量
alpha = 1;
beta  = 2;
rho   = 0.5;

% 这里直接用PPT的题目来做，这样这样对不对自己心里比较有数
% https://blog.csdn.net/peiwang245/article/details/78072130
% 大连理工大学谷俊峰
% 直接用别人的邻接矩阵
D = [ 0, 3, 1, 2;
	  3, 0, 5, 4; 
	  1, 5, 0, 2;
	  2, 4, 2, 0; ];
  
% 计算初始信息素
% 贪心算法算最短路径作为参考
% 从A出发，回到A，遍历所有点，一共是4个点，所以可以选择遍历4次，本来是5次，那么这里因为不能重复，所以最后一次是固定的
% 准备工作
k1 = 1;				% 开始从A点出发
k2 = 1;
Cnn = 0;
% 原来为0的改成无穷，这样贪心算法好算，因为不能找到自身位置
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

% 初始化信息素矩阵
tau = ones(locationNum) * (m / Cnn);			% 信息素矩阵对角线为0，其余值为m/Cnn
tau = tau - eye(locationNum) * (m / Cnn);		% 这里用全1矩阵减去对角矩阵得到




