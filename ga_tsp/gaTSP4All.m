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

%% �Ŵ��㷨����
% ������ʼ�����������
N = 50;
L = locationNum;

pc = 0.8;							% �������
pm = 0.25;							% �������

iter = 200;					% ��������
for i = 1 : N
	dna(i, :) = randperm(10);								% �������ظ����������
end

for gen = 1 : iter
	
	% ����
	% ��Ϊ�����ȫ0�����飬��߼Ӹ�����
	x1 = dna;
	x2 = dna;
	for i = 1: N
		if rand < pc							
			d = randi(N);						% ȷ����һ������ĸ���
			m = dna(d,:);						% ȷ����һ������ĸ���
			d = randi(L-1);						% ȷ������ϵ�
			x1(i,:) = [dna(i,1:d), m(d+1:L)];	% �¸��� 1������ע��һ�£���Զ��ԭʼ��Ⱥ��i������͵���������廻        
			x2(i,:) = [m(1:d), dna(i,d+1:L)];	% �¸��� 2
		end
	end
	
	% ����
	x3 = dna;
	for i = 1: N								% �����������ߵı��첻��rand�ˣ���Ϊ��Ų����ظ���������������λ���໥����
		if rand < pm
			k1 = randi([1, L]);
			k2 = randi([1, L]);
			t = x3(i, k1);
			x3(i, k1) = x3(i, k2);
			x3(i, k2) = t;
		end
	end
	
	
	% �������̭
	dna = [dna; x1; x2; x3];							% �ϲ��¾ɻ���
	
	% ������Ӧ��
	% ��ʵ����û��ô���ӣ����ǲ��
	dna = [dna, dna(:, 1) ];							% ������ص���ʼ�ĵط�
	for i = 1 : N * 4
		fi(i, 1) = 0;
		for j = 2 : L + 1
			k1 = dna(i, j);
			k2 = dna(i, j - 1);
			fi(i, 1) = fi(i, 1) + costmap(k1, k2);
		end
	end
	
	dna = [dna, fi];
	dna = sortrows(dna, L + 2);							% ����Ӧ�Ƚ���������ע����������Ǵӵ͵���
	
	% ��Ȼѡ��
	% ����Ӵ�����̫���ˣ��Ϳ�ʼ��̭
	while size(dna, 1) > N
		d = randi(size(dna, 1));
		if rand < (d - 1) / size(dna, 1)				% ��������Խ�������̭����Խ��
			dna(d,:) = [];
			fi(d, :) = [];
		end
	end
	dna = dna(:, 1:L);
	
end

route = [ dna(1, :), dna(1, 1) ];
for i = 2 : L + 1
	posK1 = route(:, i - 1);
	posK2 = route(:, i);
	
	posX = pos(posK1, 1);
	posY = pos(posK1, 2);
	posU = pos(posK2, 1);
	posV = pos(posK2, 2);
	
	quiver(posX, posY, posU - posX, posV - posY);
	
end

%% ���
%  ���������ʵ�ĵ���󣬸Ķ��ˣ�����Ҳ�ĳ�������
% �ᷢ�� �õ������ݶ���
% [8,8,8,3,3,3,3,8,8,8]
% ������һ���ط�תȦ�ģ�������Ϊ������ʱ����Щ����û�п��ǵ�

% ����1 2 3 4 5 6��6 5 4 3 2 1����
% ���ܱ�Ϊ 1 2 3 3 2 1��6 5 4 4 5 6

