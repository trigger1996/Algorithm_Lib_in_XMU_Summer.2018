clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, xȡֵ��Χ[-1, 2]

interval_max = 2;
interval_min = -1;

iter_max = 1000;		% ��������

m = 50;					% ��������
alpha = 5;				% ��Ϣ�ص������ӣ������������ת�Ƹ��ʵ�
%beta  = 1;				% ·�̵������ӣ��������ת�Ƹ��ʵģ����������������������������
rho = 0.9;				% ��Ϣ�ػӷ�����

%dx  = (interval_max - interval_min) * 0.1;
						% �ֲ��������䣬����������汻������

%
% �����������䶯

% ת�Ƹ����У�·����Ϊ��������ʹ��
% ͬʱ����Ϣ��Խǿ�ҵĵط���ת�Ƹ���Խ�ͣ�tau�ñ������ֵ��ȥ��ǰֵ���ٹ�һ��
% ����P = (tau_max - tau) / sum(tau);
% ���淢����ô��P���С����̫��
% ����P = (tau_max - tau) / tau_max;

% ��Ϣ�صļ��㣬��������ֱ��ʹ�ú���ֵ����
% tau = (1 - rho) * tau + F(x);
% ���ϵ��ƶ�ʱ����ģ�������ת�Ƹ���С�������������ƶ���������ȫ������ƶ�

% �������ϵĳ�ʼλ��
x = rand(1, m) * (interval_max - interval_min) + interval_min;		% �� [������Сֵ, �������ֵ] ��Χ���������������Ϊrand��0~1�����Կ�����ô��
x_init = x;															% ��¼��ʼֵ��������ʼ��

% ��ʼ����Ϣ�ؾ���
tau = f(x);															% ��ʼ����Ϣ�ؾ���f(x)

for iter = 1 : iter_max
	% ����ת�Ƹ���
	tau_max = max(tau);
	P = (tau_max - tau).^alpha / tau_max;							% (tau_max - tau) / sum(tau)��������������ȡֵ̫С�ˣ����Ըɴ�������ֵ
	
	% ���Ͽ�ʼ��
	for i = 1 : m
		if rand < P(i)												% ����ת�Ƹ�����ȫ������
			pos_neg = sign(rand-0.5);								% �����ϵİ취�����������������
			x(i) = x(i) + pos_neg * (interval_max - interval_min) * rand;
		else														% ����ֲ�����
			pos_neg = sign(rand-0.5);								% �����ϵİ취�����������������
			x(i) = x(i) + pos_neg * 1 / iter * rand;				% x = x + pos_neg * dx * rand������ط���ܴ�̶�Ӱ���㷨������
																	% x = x + pos_neg * 1 / iter * rand
																	% ԭ�Ĳ��ú��ߣ������Ժܺã��������ǰ�ߣ�����ȫ������
		end
	end
	
	for i = 1 : m
		if x(i) < interval_min
			x(i) = interval_min;
		end
		if x(i) > interval_max
			x(i) = interval_max;
		end
	end
	
	% ������Ϣ��
	tau = (1 - rho) * tau +  f(x);
	
end

figure(1); 
fplot(f, [-1, 2]); hold on
plot(x_init, f(x_init),'ko','linewidth',2)								% ������ʼ���λ��

figure(2);
fplot(f, [-1, 2]); hold on
plot(x, f(x),'ro','linewidth',3)									% �������ս��λ��
