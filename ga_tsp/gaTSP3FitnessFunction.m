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

% ������ʼ�����������
N = 50;
L = locationNum;
for i = 1 : N
	dna(i, :) = randperm(10);								% �������ظ����������
end

% ������Ӧ�Ⱥ���
% ��ʵ����û��ô���ӣ����ǲ��
dna_t = [dna, dna(:, 1) ];									% ������ص���ʼ�ĵط�
for i = 1 : N
	fi(i, 1) = 0;
	for j = 2 : L + 1
		k1 = dna_t(i, j);
		k2 = dna_t(i, j - 1);
		fi(i, 1) = fi(i, 1) + costmap(k1, k2);
	end
end
