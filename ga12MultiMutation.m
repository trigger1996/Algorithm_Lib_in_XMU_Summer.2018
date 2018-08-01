% ������
clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, xȡֵ��Χ[-1, 2]
%f = @(x) sin(x) + x .* cos(x);		% ��������ʽ

pc = 0.8;							% �������
pm = 0.05;							% �������

N = 50;								% ��Ⱥ����
L = 10;								% ���򳤶�
iter = 100;							% ��������
dna = randi([0, 1], [N, L]);		% ����rand()����1��ÿ������ȡֵ��Χ������2�Ǿ��������ֵ

dcd = [ 512; 256; 128; 64; 32; 16; 8; 4; 2; 1 ];
%num = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;

% �ܵ���
for gen = 1 : iter
	% ����
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
	% ���������
	x3 = dna;
	mutation_pts = randi([1, L], [1, 3]);				% ��ʱֻ��3������죬����һ����Ϊ10�Ļ�����˵��3�����������
	mutation_pts = sort(mutation_pts);
	for i = 1: N										% �������
		if rand < pm
			for j = 1: size(mutation_pts, 2)
				x3(i, mutation_pts(j)) = randi([0, 1]);
			end
		end
	end
	
	% �������̭
	dna = [dna; x1; x2; x3];							% �ϲ��¾ɻ���
	fi = f(dna * dcd * ((2 - -1) / (2^10 - 1)) - 1);	% ������Ӧ��
	dna = [dna, fi];
	dna = flipud(sortrows(dna, L + 1));					% ����Ӧ�Ƚ�������
	
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
x = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;			% ��������Ⱥ����

fplot(f, [-1, 2]); hold on;								% ��������ͼ��
plot(x, f(x),'ro','linewidth',3)						% �������ս��λ��