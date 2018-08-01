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

% TSP问题的基因长必须和地点数相同
% 而且二进制编码很难做TSP问题，最好用浮点编码或者十进制编码
N = 50;
L = locationNum;
for i = 1 : N
	dna(i, :) = randperm(10);								% 产生不重复的随机数列
end

