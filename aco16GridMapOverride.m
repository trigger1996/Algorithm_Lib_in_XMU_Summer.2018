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

iter_max = 100;

m = 50;
alpha = 1;
beta  = 7;
rho   = 0.3;
Q = 1;

startPos = 1;
endPos   = AO_Length * AO_Length;									% 为了方便邻接矩阵表示，这里的坐标不使用二维坐标
																	% 而改用0~400

% 计算邻接矩阵
l=size(G,1); 
D=zeros(l*l,l*l);													% 邻接矩阵需要考虑点与每个点之间的联系，所以必须是原坐标值的长*宽
for i=1:l
	for j=1:l
		if G(i,j)==0												% 当前点可通行才进行计算
			for m=1:l
				for n=1:l											% 遍历所有点，查找临近点
					if G(m,n)==0									% 目标点可通行才进行计算
						im=abs(i-m);jn=abs(j-n);					% 计算dx和dy
						if im+jn==1||(im==1&&jn==1)					% x、y差值为1或者x、y差值均为1才是邻近的
							D((i-1)*l+j,(m-1)*l+n)=(im+jn)^0.5;		% 这种情况下才使用勾股定理计算运动开销
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

			% 方向约束
			flag_remove = logical(zeros(size(allowList, 2), 1));
			for j = 1 : size(allowList, 2)
				currentGrid = antTraj(i, currentPos);
				nextGrid    = allowList(:, j);
				% 单坐标反解(x, y)上的坐标
				x0 = mod(currentGrid, AO_Length);
				y0 = fix(currentGrid / AO_Length) + 1;
				x1 = mod(nextGrid, AO_Length);
				y1 = fix(nextGrid / AO_Length) + 1;
				
%				% 只允许向右或向下运动
% 				if (x1 < x0 || y1 < y0)
% 					flag_remove(j, :) = 1;							% 这边只标记，在下面再删除
% 				end
				% 改为上下不约束，只允许从左到右
				%if (x1 < x0)
				%	flag_remove(j, :) = 1;
				%end
			end
			Pt(flag_remove) = [];									% 删除那些不满足方向约束的
			allowList(flag_remove') = [];
			
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
			
% 			Pr = [ allowList', Pt ];
% 			Pr = flipud(sortrows(Pr, 2));							% 完成排名
% 			if currentPos == size(antTraj, 2)
% 				antTraj = [ antTraj, zeros(m, 1) ];					% 预分配下一列全为0，为了保证如果排名选不中，补救的if在运行的时候的稳定性
% 			end
% 			for j = 1 : size(Pr, 1)
% 				if rand < Pr(j, 2)
% 					antTraj(i, currentPos + 1)  = Pr(j, 1);					
% 				end
% 			end
% 			if antTraj(i, currentPos + 1) == 0						% 排名这么写是可能选不中的，这边加个if补救一下
% 				antTraj(i, currentPos + 1)  = Pr(j, 1);
% 			end

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

%bestTraj = antTraj(1, :);		% 读取最好的那只蚂蚁的路径，顺便把最后一个0去了
%bestTraj = bestTraj(1 : size(bestTraj, 2) - 1);

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

