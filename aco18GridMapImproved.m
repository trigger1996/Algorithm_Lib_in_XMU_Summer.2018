clear;
clc;
clear all;

G=[0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 1 1 1 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 1 1 1 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0;
   0 1 1 1 0 0 1 1 1 0 1 1 1 1 0 0 0 0 0 0; 
   0 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0; 
   0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0; 
   1 1 1 1 0 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0; 
   1 1 1 1 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 1 1 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0;];
G = logical(G);														% �ĳ��߼��͵ļ�������

AO_Length = size(G, 1);												% ���򳤶�

L = AO_Length;
startGrid = [ 0, 0 ];
endGrid   = [ L, L ];

iter_max = 200;

m = 100;
alpha = 1;
beta  = 10;
rho   = 0.3;
Q = 1;

startPos = 1;
endPos   = AO_Length * AO_Length;									% Ϊ�˷����ڽӾ����ʾ����������겻ʹ�ö�ά����
																	% ������0~400

% �����ڽӾ���
G = ~G;
G_Grad = zeros(L);
for y = 1 : L
	for x = 1 : L
		if G(y, x) == 0
			G_Grad(y, x) = 0;
		else
			dx = x - endGrid(1);
			dy = y - endGrid(2);
			G_Grad(y, x) = sqrt(dx * dx + dy * dy);
		end
	end
end
G_Grad(L, L) = 0.01;												% �����и�BUG, ������������0

l=size(G,1); 
D=zeros(l*l,l*l);													% �ڽӾ�����Ҫ���ǵ���ÿ����֮�����ϵ�����Ա�����ԭ����ֵ�ĳ�*��
for i=1:l
	for j=1:l
		if G_Grad(i,j)~=0												% ��ǰ���ͨ�вŽ��м���
			for m=1:l
				for n=1:l											% �������е㣬�����ٽ���
					if G_Grad(m,n)~=0									% Ŀ����ͨ�вŽ��м���
						im=abs(i-m);jn=abs(j-n);					% ����dx��dy
						if im+jn==1||(im==1&&jn==1)					% x��y��ֵΪ1����x��y��ֵ��Ϊ1�����ڽ���
							%D((i-1)*l+j,(m-1)*l+n)=(im+jn)^0.5;
							D((i-1)*l+j,(m-1)*l+n) = (G_Grad(i,j) + G_Grad(m,n)) / 2.;
						end
					end
				end
			end
		end
	end
end
% ������������
Eta = 1 ./ D;
Eta(find(Eta == inf | isnan(Eta))) = 0;								%  1 / 0 = inf
Eta = Eta .* 100;

G = ~G;

% �����ʼ��Ϣ��
Cnn = sqrt(AO_Length * AO_Length + AO_Length * AO_Length);			%  ̰���㷨�����̫���ˣ�����ֱ����ֱ�߾������

% ��ʼ����Ϣ�ؾ���
tau = ones(AO_Length * AO_Length) * (m / Cnn);						% ��Ϣ�ؾ���Խ���Ϊ0������ֵΪm/Cnn
tau = tau - eye(AO_Length * AO_Length) * (m / Cnn);					% ������ȫ1�����ȥ�ԽǾ���õ�

% ��������
% �ҵ���·�����к�
Traj_out = cell(1,1);												% ��1��������켣����2�������濪��
Traj_out_cxx = [];
traj_out_seq = 1;

% ��ѭ������
for iter = 1 : iter_max
	
	% ������
	antTraj = zeros(m, 1);
	for i = 1 : m
		antTraj(i, :) = startPos;
	end
	% �������ϵĳ�ʼ·��
	Cxx = zeros(m, 1);
	
	for i = 1 : m
		
		% ÿֻ��������ָ����·��
		currentPos = 1;
		
		while isempty(find(antTraj(:, currentPos) == endPos))
			rows = antTraj(i, currentPos);							% �ڽӾ����еĵ�ǰ��
			allowList = find(D(rows, :) ~= 0);						% ���ڽӾ���ǰ���������õ���һ������ǰ���ĸ��
			allowList(ismember(allowList, antTraj(i, :))) = [];		% ������Щ����п��ܰ�����һ���߹��ģ����Ա���ɾ����һ���߹���
																	
			Pt = zeros(size(allowList, 2), 1);
			for j = 1 : size(allowList, 2)							% ����ת�Ƹ���
				currentGrid = antTraj(i, currentPos);
				nextGrid    = allowList(:, j);

				Pt(j, :) = tau(currentGrid, nextGrid)^alpha + Eta(currentGrid, nextGrid)^beta;
			end
			
			Pt = Pt ./ sum(Pt);
			
			if isempty(allowList)									% �����ǻ��߽�����ͬ�ģ������ǵ�����µķ������ϻ��������ͬ��Ȼ��ĵ�
				break;
			end

			% ���̶�����һ���ƶ��Ĺ켣
			% �������е�ë����������̫ǿ��
			Pr = cumsum(Pt);
			Pr = Pr ./ Pr(end);
			selectSeq = find(Pr>=rand); 
			selectSeq = min(selectSeq);
			antTraj(i, currentPos + 1) = allowList(selectSeq);
			
			currentPos = currentPos + 1;
		end
		
		% ����ÿֻ���ϵ�����
		antTraj(i, 1 : size(antTraj, 2) + 1) = [ antTraj(i, :), 0 ];% ��0��Ϊ�˽��м���
		row_end = min(find(antTraj(i, :) == 0));
		for j = 2 : row_end - 1
			currentGrid = antTraj(i, j - 1);
			nextGrid    = antTraj(i, j);
			Cxx(i, :) = Cxx(i, :) + D(currentGrid, nextGrid);
		end
	end
		
	% ������Ϣ��
	tau = (1 - rho) .* tau;											% ����
	antTraj(:, find(sum(antTraj, 1) == 0)) = [];					% ����ȫ0��
	antTraj = [antTraj, zeros(m, 1) ];								% ��һ��0�������ҵ�ÿһ�е�ĩβ
	for i = 1 : m
		row_end = min(find(antTraj(i, :) == 0));
		
		if antTraj(i, row_end - 1) ~= endPos
			continue;
		end
		for j = 2 : row_end - 1
			currentGrid = antTraj(i, j - 1);
			nextGrid    = antTraj(i, j);
			if Cxx(i) == 0
				Cxx(i) = 1e6;
			end
			tau(currentGrid, nextGrid) = tau(currentGrid, nextGrid) + Q / Cxx(i);
		end
	end

	% ��¼·��
	for i = 1 : m
		row_end = min(find(antTraj(i, :) == 0));
	
		if antTraj(i, row_end - 1) == endPos
			Traj_out{traj_out_seq, 1} = antTraj(i, 1 : row_end - 1);
			Traj_out_cxx(traj_out_seq, :) = Cxx(i);
			traj_out_seq = traj_out_seq + 1;
		end
		
	end
		
end

% ��ͼ
a = 1;																			%С�������صı߳�

% ��������
figure(1)
MM = AO_Length;
axis([0,MM,0,MM]) 
for i=1:MM 
	for j=1:MM 
		if G(i,j)==1 
			x1=j-1;y1=MM-i; 
			x2=j;y2=MM-i; 
			x3=j;y3=MM-i+1; 
			x4=j-1;y4=MM-i+1; 
			fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]); 
			hold on 
		else 
			x1=j-1;y1=MM-i; 
			x2=j;y2=MM-i; 
			x3=j;y3=MM-i+1; 
			x4=j-1;y4=MM-i+1; 
			fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]); 
			hold on 
		end 
	end 
end

% ������õ���ֻ���ϵĹ켣
index = find(Traj_out_cxx == min(Traj_out_cxx));
bestTraj = Traj_out{index, 1};

x = zeros(1, size(bestTraj, 2));
y = zeros(1, size(bestTraj, 2));
for i=1:size(bestTraj, 2) 
	x(i)=a*(mod(bestTraj(i),MM)-0.5); 
	if x(i)==-0.5 
		x(i)=MM-0.5; 
	end 
		y(i)=a*(MM+0.5-ceil(bestTraj(i)/MM)); 
end 
plot(x,y) 


% ������������
%figure(2);

