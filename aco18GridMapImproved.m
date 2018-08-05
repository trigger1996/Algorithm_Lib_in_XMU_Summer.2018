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
G = logical(G);														% 改成逻辑型的加速运算

AO_Length = size(G, 1);												% 区域长度

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
endPos   = AO_Length * AO_Length;									% 为了方便邻接矩阵表示，这里的坐标不使用二维坐标
																	% 而改用0~400

% 计算邻接矩阵
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
G_Grad(L, L) = 0.01;												% 这里有个BUG, 这个点算出来是0

l=size(G,1); 
D=zeros(l*l,l*l);													% 邻接矩阵需要考虑点与每个点之间的联系，所以必须是原坐标值的长*宽
for i=1:l
	for j=1:l
		if G_Grad(i,j)~=0												% 当前点可通行才进行计算
			for m=1:l
				for n=1:l											% 遍历所有点，查找临近点
					if G_Grad(m,n)~=0									% 目标点可通行才进行计算
						im=abs(i-m);jn=abs(j-n);					% 计算dx和dy
						if im+jn==1||(im==1&&jn==1)					% x、y差值为1或者x、y差值均为1才是邻近的
							%D((i-1)*l+j,(m-1)*l+n)=(im+jn)^0.5;
							D((i-1)*l+j,(m-1)*l+n) = (G_Grad(i,j) + G_Grad(m,n)) / 2.;
						end
					end
				end
			end
		end
	end
end
% 计算启发因子
Eta = 1 ./ D;
Eta(find(Eta == inf | isnan(Eta))) = 0;								%  1 / 0 = inf
Eta = Eta .* 100;

G = ~G;

% 计算初始信息素
Cnn = sqrt(AO_Length * AO_Length + AO_Length * AO_Length);			%  贪心算法算这个太大了，所以直接用直线距离代替

% 初始化信息素矩阵
tau = ones(AO_Length * AO_Length) * (m / Cnn);						% 信息素矩阵对角线为0，其余值为m/Cnn
tau = tau - eye(AO_Length * AO_Length) * (m / Cnn);					% 这里用全1矩阵减去对角矩阵得到

% 其他变量
% 找到出路的序列号
Traj_out = cell(1,1);												% 第1列用来存轨迹，第2列用来存开销
Traj_out_cxx = [];
traj_out_seq = 1;

% 主循环迭代
for iter = 1 : iter_max
	
	% 放蚂蚁
	antTraj = zeros(m, 1);
	for i = 1 : m
		antTraj(i, :) = startPos;
	end
	% 清零蚂蚁的初始路径
	Cxx = zeros(m, 1);
	
	for i = 1 : m
		
		% 每只蚂蚁运行指定的路程
		currentPos = 1;
		
		while isempty(find(antTraj(:, currentPos) == endPos))
			rows = antTraj(i, currentPos);							% 邻接矩阵中的当前行
			allowList = find(D(rows, :) ~= 0);						% 在邻接矩阵当前行搜索，得到下一个可以前往的格点
			allowList(ismember(allowList, antTraj(i, :))) = [];		% 但是这些格点中可能包括上一次走过的，所以必须删掉上一次走过的
																	
			Pt = zeros(size(allowList, 2), 1);
			for j = 1 : size(allowList, 2)							% 计算转移概率
				currentGrid = antTraj(i, currentPos);
				nextGrid    = allowList(:, j);

				Pt(j, :) = tau(currentGrid, nextGrid)^alpha + Eta(currentGrid, nextGrid)^beta;
			end
			
			Pt = Pt ./ sum(Pt);
			
			if isempty(allowList)									% 蚂蚁是会走进死胡同的，这行是调完底下的发现蚂蚁会钻进死胡同，然后改的
				break;
			end

			% 轮盘赌算下一次移动的轨迹
			% 排名法有点毛病，收敛性太强了
			Pr = cumsum(Pt);
			Pr = Pr ./ Pr(end);
			selectSeq = find(Pr>=rand); 
			selectSeq = min(selectSeq);
			antTraj(i, currentPos + 1) = allowList(selectSeq);
			
			currentPos = currentPos + 1;
		end
		
		% 计算每只蚂蚁的消耗
		antTraj(i, 1 : size(antTraj, 2) + 1) = [ antTraj(i, :), 0 ];% 补0，为了进行计算
		row_end = min(find(antTraj(i, :) == 0));
		for j = 2 : row_end - 1
			currentGrid = antTraj(i, j - 1);
			nextGrid    = antTraj(i, j);
			Cxx(i, :) = Cxx(i, :) + D(currentGrid, nextGrid);
		end
	end
		
	% 更新信息素
	tau = (1 - rho) .* tau;											% 蒸发
	antTraj(:, find(sum(antTraj, 1) == 0)) = [];					% 清零全0列
	antTraj = [antTraj, zeros(m, 1) ];								% 加一列0，用来找到每一行的末尾
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

	% 记录路径
	for i = 1 : m
		row_end = min(find(antTraj(i, :) == 0));
	
		if antTraj(i, row_end - 1) == endPos
			Traj_out{traj_out_seq, 1} = antTraj(i, 1 : row_end - 1);
			Traj_out_cxx(traj_out_seq, :) = Cxx(i);
			traj_out_seq = traj_out_seq + 1;
		end
		
	end
		
end

% 作图
a = 1;																			%小方格象素的边长

% 环境构造
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

% 作出最好的那只蚂蚁的轨迹
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


% 作出收敛曲线
%figure(2);

