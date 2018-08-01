clear;
clc;
clear all;

%
% 构建整个随机地图，然后通过地图来计算目标值
%
locationNum = 10;
pos = randi([0, 10000], [locationNum, 2]) / 100.;			% 产生一个0~100的(x, y)坐标向量，randi的参数前者是范围，是矩阵的Size

figure;
for i = 1 : locationNum
	plot(pos(i, 1), pos(i, 2),'ro','linewidth',3); hold on;						% 画出最终解的位置
end

costmap = zeros(10, 10);
for i = 1 : locationNum
	for j = 1 : locationNum
		dx = pos(i, 1) - pos(j, 1);
		dy = pos(i, 2) - pos(j, 2);
		costmap(i, j) = sqrt(dx*dx + dy*dy);
	end
end

