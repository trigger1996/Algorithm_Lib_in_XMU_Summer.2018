clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, xȡֵ��Χ[-1, 2]
%f = @(x) sin(x) + x .* cos(x);		% �������ʽ

pc = 0.8;							% �������
pm = 0.05;							% �������

N = 50;								% ��Ⱥ����
L = 10;								% ���򳤶�
iter = 6;							% ��������
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
	x3 = dna;
	for i = 1: N								% �������
		if rand < pm
			x3(i,randi(L)) = randi([0, 1]);
		end
	end
	
	% �������̭
	% ��������ѡ��
	% ������������һ�����ѡ�����������б�����Ȼ����Ӧ����ߵ�����ȥ
	% ���������ظ�ѡ��
	% д����������
	dna_t = [dna; x1; x2; x3];							% �ϲ��¾ɻ���
	fi = f(dna_t * dcd * ((2 - -1) / (2^10 - 1)) - 1);	% ������Ӧ��
	dna_t = [dna_t, fi];
	
	nIden = size(dna_t, 1);
	nSel  = N;
	dna   = [];
	for i = 1 : nSel
		tournSel = randi([1, nIden], [1, 5]);				% ÿ�������ѡ�������ٸ������о���
		tournSel = sort(tournSel);							% ������Ӧ�Ⱦ��Ǳ�����
		
		compti = dna_t(tournSel, :);						% ѡ�е��������
		compti = flipud(sortrows(compti, L + 1));			% ����Ӧ�Ƚ�������
		dna(i, :) = compti(1, :);							% ��õ�һ����ѡ��
		
	end
	dna = dna(:, 1:L);	
end
x = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;			% ��������Ⱥ����

fplot(f, [-1, 2]); hold on;								% ��������ͼ��
plot(x, f(x),'ro','linewidth',3)						% �������ս��λ��
