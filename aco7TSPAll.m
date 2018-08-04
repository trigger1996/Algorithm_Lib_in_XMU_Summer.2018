clear;
clc;
clear all;

%%
%  ��һ���ֵĹ���

iter_max = 100;

locationNum = 4;
cityList    = 1 : locationNum;

m = 6;			% ��������, 4
alpha = 1;
beta  = 2;
rho   = 0.5;

% ����ֱ����PPT����Ŀ���������������Բ����Լ�����Ƚ�����
% https://blog.csdn.net/peiwang245/article/details/78072130
D = [ 0, 3, 1, 2;
	  3, 0, 5, 4; 
	  1, 5, 0, 2;
	  2, 4, 2, 0; ];
  
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
Cnn = Cnn +D(k1, 1);							% ���Ҫ�ӵ�ǰ��ص�A�㣬��Ϊ���Ƶľ����ƻ��ò���ˣ�����ֱ����ԭ������

D0 = D;											% ǿ��֢

% ������Ҫע�⣬���������ӵ�ʱ����Ҫ������������߱��봦������е�0
D(find(D == 0)) = 1e-6;
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
			% ���̶�����һ���ƶ��ĳ���
			Pr = cumsum(Pt);											% ���ۼƸ���
			seq = find(rand <= Pr);
			seq = min(seq);												% ��������ָ�룬����ָ�����պñ��ۼƸ���С���ۼƸ����Դ�һ��ľ���ָ��λ��
			antTraj(j, i + 1)  = allowList(seq);						% �켣д��
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

end

% ��������ֻ·����̵����ϵ�·������
seq_out = find(min(Cxx));
bestTraj = antTraj(seq_out, :);

