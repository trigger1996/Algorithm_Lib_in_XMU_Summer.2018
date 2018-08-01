clear;
clc;
clear all;

%
% �ڵ�ͼ������ϵĻ�����
% ���е�ͼ����
%
locationNum = 10;
pos = randi([0, 10000], [locationNum, 2]) / 100.;			% ����һ��0~100��(x, y)����������randi�Ĳ���ǰ���Ƿ�Χ���Ǿ����Size

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
	for i = 1: N
		d1 = randi(N);
		d2 = randi(N);
		% ����Ϊ��͵����ֱ��ֲ���˱�ĺ�����swap�����intercross.m
		a = dna(d1, :);
		b = dna(d2, :);
		if rand < pc							
			r1=randsrc(1,1,[1:L]);
			r2=randsrc(1,1,[1:L]);
			if r1~=r2
				a0=a;b0=b;							% �ܵĽ���֮ǰ����һ��
				s=min([r1,r2]);
				e=max([r1,r2]);
				for i=s:e
					a1=a;b1=b;						% ÿһ�ֽ���֮ǰ����һ��
					a(i)=b0(i);						% ��������
					b(i)=a0(i);
					x=find(a==a(i));				% ���������꿴һ�£���ǰ�����У���û�кͽ�������ظ���
					y=find(b==b(i));
					i1=x(x~=i);						% ��Ȼ��find�ҵ��Լ������ǳ����Լ���������б����ͬ�ģ�����i1��¼������
					i2=y(y~=i);						% ����ÿ��ֻ����һ������������2����ͬ�ģ�һ�����Լ���һ����ԭ�еģ��������ֻ���Լ�����ôi1Ϊ�գ�ͬ����i2Ҳ������
					if ~isempty(i1)					% ���i1�ǿ�
						a(i1)=a1(i);				% һģһ����������������һ���Ǹ��滻�ģ�һ����ԭ�����еģ���ô��ȡ�޸��Ǹ�ԭ�����еģ�ԭ�����еĸĳ�ǰ�ߣ������滻��֮ǰ�����֣������γɽ���
					end
					if ~isempty(i2)
						b(i2)=b1(i);
					end
				end
			end
		end
		x1(i, :) = a;
		x2(i, :) = b;
	end
	
	% ����
	x3 = dna;
	for i = 1: N								% �����������ߵı��첻��rand�ˣ���Ϊ��Ų����ظ���������������λ���໥����
		if rand < pm
			k1 = randi([1, L]);
			k2 = randi([1, L]);
			t = x3(i, k1);
			x3(i, k1) = x3(1, k2);
			x3(i, k2) = t;
		end
	end
	
	
	% �������̭
	dna_t = [dna; x1; x2; x3];							% �ϲ��¾ɻ���
	
	% ������Ӧ��
	% ��ʵ����û��ô���ӣ����ǲ��
	dna_t = [dna, dna(:, 1) ];							% ������ص���ʼ�ĵط�
	for i = 1 : N
		fi(i, 1) = 0;
		for j = 2 : L + 1
			k1 = dna_t(i, j);
			k2 = dna_t(i, j - 1);
			fi(i, 1) = fi(i, 1) + costmap(k1, k2);
		end
	end
	
	dna = [dna, fi];
	dna = sortrows(dna, L + 1);							% ����Ӧ�Ƚ���������ע����������Ǵӵ͵���
	
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
	
	% ��¼fi������ֵ�Բ鿴�������
	fi_best(gen) = fi(1, :);
	
end

%% ��ͼ

% ����λ��
figure;
for i = 1 : locationNum
	plot(pos(i, 1), pos(i, 2),'ro','linewidth',3); hold on;	% �������ս��λ��
	text(pos(i, 1) + 1, pos(i, 2) + 1, num2str(i), 'FontSize', 12);
end

% ����·��
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

% �����������
figure;
plot(1 : 200, fi_best);
