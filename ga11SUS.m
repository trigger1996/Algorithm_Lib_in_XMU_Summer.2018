clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, xȡֵ��Χ[-1, 2]
%f = @(x) sin(x) + x .* cos(x);		% �������ʽ

pc = 0.8;							% �������
pm = 0.1;							% �������

N = 50;								% ��Ⱥ����
L = 20;								% ���򳤶�
iter = 150;							% ��������
dna = randi([0, 1], [N, L]);		% ����rand()����1��ÿ������ȡֵ��Χ������2�Ǿ��������ֵ

dcd = [ 524288; 262144; 131072; 65536; 32768; 16384; 8192; 4096; 2048; 1024; 512; 256; 128; 64; 32; 16; 8; 4; 2; 1 ];
%num = dna * dcd * ((2 - -1) / (2^10 - 1)) - 1;

% �ܵ���
for gen = 1 : iter
	% ����
	for i = 1 : N
		if rand < pc
			d  = randi(N);							% ��ȡ1�����ڽ����ĸ���
			m  = dna(d, :);

			d1 = randi(L);
			d2 = randi(L);							% ����õ������ϵ�
			if d1 > d2
				t = d1;								% ����������������ر�֤�ϵ�d1��d2֮ǰ
				d1 = d2;
				d2 = t;
				clear t; 
			end
			% ��������ֹ������
			if d1 >= 20
				d1 = 18;
				d2 = 19;
			end
			if d2 >= 20
				d2 = 19;
			end

			x1(i, :) = [ dna(i, 1 : d1-1), m(d1 : d2),       dna(i, d2+1 : L) ];
			x2(i, :) = [ m(1 : d1-1),      dna(i, d1 : d2),  m(d2+1 : L) ];
		end
	end
	
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
	fi = f(dna_t * dcd * ((2 - -1) / (2^20 - 1)) - 1);	% ������Ӧ��
	%dna = [dna, fi];

	% ������Ȳ�����SUS
	% ��һ������ֵ
	fi_min = min(fi);
	fi = fi - fi_min;
	fi = cumsum(fi);									% �ۼ�
	fi = fi / fi(size(fi, 1));							% �ۼƸ������һ��ֵһ�������ֵ����ȫ����ӣ����������͹�һ����
	
	% ����ָ��
	% ָ���Ǿ��ȸ����ģ����Լ�����Ȼ�ǣ�
	% ��ǰ��Ⱥ�� / Ŀ����Ⱥ�����������ٸ�ȡһ��
	% ���Ƕ��ڸ��ʣ����ڸ�����0~1���������ֵҲ����������϶�Ҫ��һ����������Ϊ����1��ȡ���ٸ�������
	% ����� 1 / nSel
	nInd = size(fi, 1);
	nSel = N;
	interval = 1 / nSel;
	susPtr = 0. : interval : 1.;
	susPtr = susPtr + rand / nSel;		% rand / nSel < 1 / nSel��rand / nSel��һ��С�����䲽���������
	
	% ����ָ��ȡ��
	j = 1;
	dna = [];
	for i = 1 : nInd
		if j > nSel
			break;
		end
		if fi(i) > susPtr(j)
			dna(j, :) = dna_t(i, :);
			j = j + 1;
		end
	end
	
end
x = dna * dcd * ((2 - -1) / (2^20 - 1)) - 1;			% ��������Ⱥ����

fplot(f, [-1, 2]); hold on;								% ��������ͼ��
plot(x, f(x),'ro','linewidth',3)						% �������ս��λ��
