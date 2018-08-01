clear;
clc;
clear all;

%
% �������������ͼ��Ȼ��ͨ����ͼ������Ŀ��ֵ
%
locationNum = 10;
pos = randi([0, 10000], [locationNum, 2]) / 100.;			% ����һ��0~100��(x, y)����������randi�Ĳ���ǰ���Ƿ�Χ���Ǿ����Size

figure;
for i = 1 : locationNum
	plot(pos(i, 1), pos(i, 2),'ro','linewidth',3); hold on;						% �������ս��λ��
end

costmap = zeros(10, 10);
for i = 1 : locationNum
	for j = 1 : locationNum
		dx = pos(i, 1) - pos(j, 1);
		dy = pos(i, 2) - pos(j, 2);
		costmap(i, j) = sqrt(dx*dx + dy*dy);
	end
end

