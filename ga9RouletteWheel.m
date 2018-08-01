clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, xȡֵ��Χ[-1, 2]
%f = @(x) sin(x) + x .* cos(x);		% �������ʽ

pc = 0.9;							% �������
pm = 0.1;							% �������

N = 125;							% ��Ⱥ����
L = 10;								% ���򳤶�
iter = 150;							% ��������
dna = randi([0, 1], [N, L]);		% ����rand()����1��ÿ������ȡֵ��Χ������2�Ǿ��������ֵ

dcd = [ 512; 256; 128; 64; 32; 16; 8; 4; 2; 1 ];
%num = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;

% �ܵ���
for gen = 1 : iter
	% ����
	for i = 1: size(dna, 1)
		if rand < pc							
			d = randi(size(dna, 1));			% ȷ����һ������ĸ���
			m = dna(d,:);						% ȷ����һ������ĸ���
			d = randi(L-1);						% ȷ������ϵ�
			x1(i,:) = [dna(i,1:d), m(d+1:L)];	% �¸��� 1������ע��һ�£���Զ��ԭʼ��Ⱥ��i������͵���������廻        
			x2(i,:) = [m(1:d), dna(i,d+1:L)];	% �¸��� 2
		end
	end
	
	% ����
	x3 = dna;
	for i = 1: size(dna, 1)						% �������
		if rand < pm
			x3(i,randi(L)) = randi([0, 1]);
		end
	end
	
	% �������̭
	% ���̶ķ�
	dna_t = [dna; x1; x2; x3];							% �ϲ��¾ɻ���
	fi = f(dna_t * dcd * ((2 - -1) / (2^10 - 1)) - 1);	% ������Ӧ��
	%dna = [dna, fi];

	% ���̶ķ�
	% ��һ������ֵ
	norm = 0;
	for i = 1 : size(dna, 1);
		norm = norm + fi(i);
	end
	fi  = fi / norm;
	% �����ۼƸ���
	pfi(1) = fi(1);
	for i = 2 : size(dna, 1);
		pfi(i) = pfi(i - 1) + fi(i);
	end
	% ת����
	dna = [];
	for i = 1 : N
		t = rand;
		for j = 1 : size(dna_t, 1)
			if t <= pfi(j)								% t < pfi(j)
				dna(i, :) = dna_t(j, :);
				break;
			end
		end
	end
	
end
x = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;			% ��������Ⱥ����

fplot(f, [-1, 2]); hold on;								% ��������ͼ��
plot(x, f(x),'ro','linewidth',3)						% �������ս��λ��
