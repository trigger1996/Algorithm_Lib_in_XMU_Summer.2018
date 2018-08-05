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

a = 1;																%С�������صı߳�


AO_Length = size(G, 1);												% ���򳤶�

iter_max = 100;

m = 50;
alpha = 1;
beta  = 7;
rho   = 0.3;
Q = 1;

startPos = 1;
endPos   = AO_Length * AO_Length;									% Ϊ�˷����ڽӾ����ʾ����������겻ʹ�ö�ά����
																	% ������0~400

% �����ڽӾ���
l=size(G,1); 
D=zeros(l*l,l*l);													% �ڽӾ�����Ҫ���ǵ���ÿ����֮�����ϵ�����Ա�����ԭ����ֵ�ĳ�*��
for i=1:l
	for j=1:l
		if G(i,j)==0												% ��ǰ���ͨ�вŽ��м���
			for m=1:l
				for n=1:l											% �������е㣬�����ٽ���
					if G(m,n)==0									% Ŀ����ͨ�вŽ��м���
						im=abs(i-m);jn=abs(j-n);					% ����dx��dy
						if im+jn==1||(im==1&&jn==1)					% x��y��ֵΪ1����x��y��ֵ��Ϊ1�����ڽ���
							D((i-1)*l+j,(m-1)*l+n)=(im+jn)^0.5;		% ��������²�ʹ�ù��ɶ�������˶�����
						end
					end
				end
			end
		end
	end
end
% ������������
Ex = a * (mod(endPos, AO_Length) - 0.5);							%��ֹ�������
if Ex == -0.5 
	Ex = AO_Length - 0.5; 
end
Ey = a * (AO_Length + 0.5 - ceil(endPos / AO_Length));				% ��ֹ��������
Eta = zeros(AO_Length * AO_Length, 1);								% ����ʽ��Ϣ��ȡΪ��Ŀ����ֱ�߾���ĵ���
%��������ʽ��Ϣ����
for i = 1 : AO_Length * AO_Length
	ix = a * (mod(i, AO_Length)-0.5); 
	if ix == -0.5 
		ix = AO_Length - 0.5; 
	end 
	iy = a * (AO_Length + 0.5 - ceil(i / AO_Length));  
	if i ~= endPos 
		Eta(i) = 1 / ((ix - Ex)^2 + (iy - Ey)^2)^0.5; 
	else 
		Eta(i) = 100;
	end 
end 


% ��ʼ����Ϣ�ؾ���
tau = ones(AO_Length * AO_Length) * 8;						% ��Ϣ�ؾ���Խ���Ϊ0������ֵΪm/Cnn
tau = tau - eye(AO_Length * AO_Length) * 8;					% ������ȫ1�����ȥ�ԽǾ���õ�

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

				Pt(j, :) = tau(currentGrid, nextGrid)^alpha + Eta(nextGrid)^beta;
			end
			
			% ����Լ��
			flag_remove = logical(zeros(size(allowList, 2), 1));
			for j = 1 : size(allowList, 2)
				currentGrid = antTraj(i, currentPos);
				nextGrid    = allowList(:, j);
				% �����귴��(x, y)�ϵ�����
				x0 = a * (mod(currentGrid, AO_Length) - 0.5); 
				if x0 == -0.5 
					x0 = AO_Length - 0.5; 
				end 
				y0 = a * (AO_Length + 0.5 - ceil(currentGrid / AO_Length)); 				

				x1 = a * (mod(nextGrid, AO_Length) - 0.5); 
				if x1 == -0.5 
					x1 = AO_Length - 0.5; 
				end 
				y1 = a * (AO_Length + 0.5 - ceil(nextGrid / AO_Length));  
		

				% ֻ���������
				if (x1 < x0)
					flag_remove(j, :) = 1;
				end
			end
			Pt(flag_remove) = [];									% ɾ����Щ�����㷽��Լ����
			allowList(flag_remove') = [];
			
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

