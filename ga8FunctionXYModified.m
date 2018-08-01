clear;
clc;
clear all;

f = @(x, y) x .* cos(2 * pi * y) + y .* sin(2 * pi * x);	% x����[-2, 2] y����[-5, 3]�����Գ�����

pc = 0.8;							% �������
pm = 0.05;							% �������

N = 200;							% ��Ⱥ����
L = 21;								% ���򳤶�
iter = 15;							% �����������޸ĵ��������Բ鿴�������
dna = randi([0, 1], [N, L]);		% ����

dcd_x = [ 512; 256; 128; 64; 32; 16; 8; 4; 2; 1; 0;    0;   0;   0;   0;  0;  0;  0; 0; 0; 0 ];
dcd_y = [ 0;   0;   0;   0;  0;  0;  0; 0; 0; 0; 1024; 512; 256; 128; 64; 32; 16; 8; 4; 2; 1 ];

for gen = 1 : iter;
	
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
			end

			x1(i, :) = [ dna(i, 1 : d1-1), m(d1 : d2),       dna(i, d2+1 : L) ];
			x2(i, :) = [ m(1 : d1-1),      dna(i, d1 : d2),  m(d2+1 : L) ];
		end
	end
	
	% ����
	x3 = dna;
	for i = 1: N
		if rand < pm
			x3(i,randi(L)) = randi([0, 1]);
		end
	end
	
	% ����
	dna = [dna; x1; x2; x3];							% �ϲ��¾ɻ���
	x = dna * dcd_x * (2 + 2) / 1023 - 2;
	y = dna * dcd_y * (3 + 5) / 2047 - 5;
	fi = f(x, y);										% ������Ӧ��
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

x = dna * dcd_x * (2 + 2) / 1023 - 2;
y = dna * dcd_y * (3 + 5) / 2047 - 5;
ezmesh(f, [-2, 2, -5, 1]); hold on;						% ��������ͼ��
plot3(x, y, f(x, y),'ro','linewidth',3)					% �������ս��λ��
