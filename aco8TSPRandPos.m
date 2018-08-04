clear;
clc;
clear all;

%%
%  ��һ���ֵĹ���

iter_max = 100;

locationNum = 8;
cityList    = 1 : locationNum;

m = ceil(locationNum * 1.5);								% ����������ȡ����������1.5������
alpha = 1;
beta  = 2;
rho   = 0.3;

% �������ĳ��������
pos = randi([0, 10000], [locationNum, 2]) / 1000.;			% ����һ��0~1000��(x, y)����������randi�Ĳ���ǰ���Ƿ�Χ���Ǿ����Size
% �����ڽӾ���
D = zeros(locationNum);
for i = 1 : locationNum
	for j = 1 : locationNum
		dx = pos(i, 1) - pos(j, 1);
		dy = pos(i, 2) - pos(j, 2);
		D(i, j) = sqrt(dx*dx + dy*dy);
	end
end

  
% �����ʼ��Ϣ��
% ̰���㷨�����·����Ϊ�ο�
k1 = 1;				% ��ʼ��A�����
k2 = 1;
Cnn = 0;

D0 = D;
D0(find(D == 0)) = inf;
% ����̰���㷨����ʵ��һ�����˹�����Ҳ�ǿ��Եģ��Ͼ������ص�
for i = 1 : size(D0, 1) - 1
	Cnn = Cnn + min(D0(k1, :));
	k2 = find(D0(k1, :) == min(D0(k1, :)));		% �ҵ�������Сֵ��Ԫ�ص��±ꡪ��Ŀ�����ҵ�̰���㷨���������һ��λ�ã���һ��λ��Ҫ�����￪ʼ����
	D0(:, k1) = inf;							% �ϻ��߹��ĵ㲻��ȥ�ˣ���������������޴�ֱ�ӽ�ֹ
	k1 = k2;									% ���ݲ�����׼����һ�ּ���
end
Cnn = Cnn + D(k1, 1);							% ���Ҫ�ӵ�ǰ��ص�A�㣬��Ϊ���Ƶľ����ƻ��ò���ˣ�����ֱ����ԭ������

D0 = D;											% ǿ��֢

% ������Ҫע�⣬���������ӵ�ʱ����Ҫ������������߱��봦������е�0
D(find(D == 0)) = 1e10;
Eta = 1 ./ D;									% �������ӣ�ȡ����ĵ������������������ת�Ƹ��ʵĵ�

% ��ʼ����Ϣ�ؾ���
tau = ones(locationNum) * (m / Cnn);			% ��Ϣ�ؾ���Խ���Ϊ0������ֵΪm/Cnn
tau = tau - eye(locationNum) * (m / Cnn);		% ������ȫ1�����ȥ�ԽǾ���õ�

%%

for iter = 1 : iter_max

	antTraj = zeros(m, locationNum);									% ���Ϲ켣
	antTraj(:, 1) = randi([1, locationNum], [1, m]);					% �������ϳ�ʼλ�ã����Ͽ�ʼʱ�ǿ��Է���ͬһ�����е�

	for i = 1 :locationNum - 1											% �ִΣ�����ԭ�ƻ��⻷�����ϵģ�Ȼ��������̫���ˣ����������ִ�
		for j = 1 : m
			allowList = cityList(~ismember(cityList, antTraj(j, :)));	% �����ǲ�����Щ������ȥ��

			Pt = [];													% �Ӹ�����
			for k = 1 : size(allowList, 2)
				currentLoc = antTraj(j, i);
				nextLoc    = allowList(k);
				Pt(k, :) = tau(currentLoc, nextLoc)^alpha + Eta(currentLoc, nextLoc)^beta;
																		% ���ݹ�ʽ����ת�Ƹ��ʣ�Pt��˳���allowList��һ����
			end
			Pt = Pt ./ sum(Pt);											% ���ݹ�ʽ�������
			% ����������ѡ����һ�ε�·��
			% ԭ�����������̶ģ�Ȼ��һֱ��������������һ������
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
	antTraj = [ antTraj, antTraj(:, 1) ];								% �켣���ϳ�ʼ�㣬���ϱ���ص���ʼ��



	% �����ڽӾ��������·��
	Cxx = zeros(m, 1);
	for i = 1 : m
		for j = 1 : locationNum
			currentLoc = antTraj(i, j);
			nextLoc    = antTraj(i, j + 1);
			Cxx(i, :) = Cxx(i, :) + D(currentLoc, nextLoc);				% �ܶ�����������������ǲ�д�ģ����Կ��������ѿ�����������ʵ���Ǹ����ڽӾ��������
		end
	end

	% ������Ϣ�ؾ���
	% ��������������������
	tau = (1 - rho) * tau;

	% Ȼ��������ϵ�·��������Ϣ�ؾ���
	% ���ݹ�ʽ��ÿ�θ��¶���ֱ�ӼӾͺ��˱��������˶�����·�̵ĵ����ͺ���
	for i = 1 : m
		for j = 1 : locationNum
			currentLoc = antTraj(i, j);
			nextLoc    = antTraj(i, j + 1);
			tau(currentLoc, nextLoc) = tau(currentLoc, nextLoc) + 1 / Cxx(i, :);
		end
	end

	% ��¼���·��
	seq_out = find(min(Cxx));
	bestTrajCxx(iter, :) = Cxx(seq_out, :);
	
end

% ��������ֻ·����̵����ϵ�·������
seq_out = find(min(Cxx));
bestTraj = antTraj(seq_out, :);

figure(1);
for i = 1 : locationNum
	plot(pos(i, 1), pos(i, 2),'ro','linewidth',3); hold on;						% ���������
	text(pos(i, 1) + .2, pos(i, 2) + .2, num2str(i), 'FontSize', 12);				% ����������
end

% ����·��
for i = 2 : locationNum + 1
	posK1 = bestTraj(:, i - 1);
	posK2 = bestTraj(:, i);
	
	posX = pos(posK1, 1);
	posY = pos(posK1, 2);
	posU = pos(posK2, 1);
	posV = pos(posK2, 2);
	
	quiver(posX, posY, posU - posX, posV - posY);
	
end

% ������������
figure(2);
plot(1 : iter, 	bestTrajCxx);

% ������·��
bestTraj
