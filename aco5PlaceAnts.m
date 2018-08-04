clear;
clc;
clear all;

%%
%  ��һ���ֵĹ���

locationNum = 4;
cityList    = 1 : locationNum;

m = 3;			% ��������
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

% ������Ҫע�⣬���������ӵ�ʱ����Ҫ������������߱��봦�������е�0
D(find(D == 0)) = 1e-6;
Eta = 1 ./ D;									% �������ӣ�ȡ����ĵ������������������ת�Ƹ��ʵĵ�

% ��ʼ����Ϣ�ؾ���
tau = ones(locationNum) * (m / Cnn);			% ��Ϣ�ؾ���Խ���Ϊ0������ֵΪm/Cnn
tau = tau - eye(locationNum) * (m / Cnn);		% ������ȫ1�����ȥ�ԽǾ���õ�

%%
% ��һ���ֵĹ������������ϣ��������

antTraj = zeros(m, locationNum);			% ���Ϲ켣
antTraj(:, 1) = randi([1, 4], [1, 3]);		% �������ϳ�ʼλ�ã����Ͽ�ʼʱ�ǿ��Է���ͬһ�����е�

for i = 1 :locationNum - 1											% �ִΣ�����ԭ�ƻ��⻷�����ϵģ�Ȼ��������̫���ˣ����������ִ�
	for j = 1 : m
		allowList = cityList(~ismember(cityList, antTraj(j, :)));	% �����ǲ�����Щ������ȥ��
																	% ismemeber�ǲ�ѯantTraj����Щ��Ա����cityList�������±�
																	% Ȼ��ȡ���Ժ������Щ���ǣ����ڹ켣�ϵľ��ǿ���ȥ��
		Pt = [];	
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