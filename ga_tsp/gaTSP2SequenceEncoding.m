clear;
clc;
clear all;

%
% �ڵ�ͼ������ϵĻ�����
% ���е�ͼ����
%
locationNum = 10;
pos = randi([0, 10000], [locationNum, 2]) / 100.;			% ����һ��0~100��(x, y)����������randi�Ĳ���ǰ���Ƿ�Χ���Ǿ����Size

figure;
for i = 1 : locationNum
	plot(pos(i, 1), pos(i, 2),'ro','linewidth',3); hold on;	% �������ս��λ��
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

% TSP����Ļ��򳤱���͵ص�����ͬ
% ���Ҷ����Ʊ��������TSP���⣬����ø���������ʮ���Ʊ���
N = 50;
L = locationNum;
for i = 1 : N
	dna(i, :) = randperm(10);								% �������ظ����������
end

