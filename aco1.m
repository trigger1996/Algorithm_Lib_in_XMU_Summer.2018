%蚂蚁算法test 
%用产生的一个圆上的十个点来检验蚂蚁算法
  
clc
clear
  
%参数
alpha = 1 ;                               %信息素指数
beta = 5  ;                                %启发指数
rho = 0.5 ;                                %挥发系数
n = 16 ;                                   %城市个数
k = 20 ;                                    %迭代次数
m = n - 1 ;                                %蚂蚁只数，这里取比城市数目少一的蚂蚁只数 n - 1
Q = 100 ;
bestr = inf ;
%产生一个圆上的十个点
x = zeros(1,n) ;
y = x ;
for i = 1 : (n/2)                         
    x(i) = rand * 20 ;
    y(i) = sqrt(100 - (x(i) - 10) ^ 2) + 10;
end
for i = (n/2 + 1) : n
    x(i) = rand * 20 ;
    y(i) = - sqrt(100 - (x(i) - 10) ^ 2) + 10;
end
%plot(x,y,'.') ;
plot(x, y,'ro','linewidth',3); hold on;
%计算距离
d = zeros(n,n) ;
for i = 1 : n
    for j = 1 : n
        d(i,j) = sqrt( ( x(i) - x(j) ) ^ 2 + ( y(i) - y(j) ) ^ 2) ;
    end
end
temp = min(d) ;
dmin = temp(1) ;
tau = ones(n,n) ;
%tau = tau ./ (n * dmin) ;                   %初始化tau信息素矩阵
  
%开始迭代
for i = 1 : k  
    %初始化
    visited = zeros(m,n) ;                  %用visited 来储存所有蚂蚁走过的城市 m×n 其中未到达的城市为0
    visited(:,1) = (randperm(n,m))';        %将m只蚂蚁随机放在n座城市 即产生一列1到n的随机数进行第一列数据的更新
    for b = 2 : n                           %所有蚂蚁都走到第b个城市时
        current = visited(:,(b-1)) ;         %所有蚂蚁现在所在城市 m×1
        allow = zeros(m,(n - b + 1)) ;
         
        for a = 1 : m
            j = 1 ;
            for s = 1 : n
                if length(find(visited(a,:) == s)) == 0
                   allow(a,j) = s ;
                   j = j + 1 ;
                end
            end
        end
         
        l = n-b+1 ;
        for a = 1 : m                       %分析第a只蚂蚁
            p = zeros(1,l) ;
            for j = 1 : l                   %根据下式来选择下一个城市
                p(j) = ( ( tau( current(a,1) , allow(a,j) ) ) ^ alpha ) * ( ( 1 / d( current(a,1) , allow(a,j) ) ) ^ beta ) ;
            end
            p = p ./ sum(p) ;               %采用轮盘赌的方式
            p = cumsum(p) ;
            pick = rand ;
            for c = 1 : l
                if pick < p(c)
                    visited(a,b) = allow(a,c) ;          %找到符合要求的城市 并 记入蚂蚁a的路径中
                    break ;
                end
            end
        end
    end
    %计算每只蚂蚁所走的路径总长
    L = zeros(1,m) ;
    for a = 1 : m
        t = d(visited(a,n),visited(a,1)) ;
        for b = 1 : (n - 1)
            t = t + d(visited(a,b),visited(a,(b + 1)));
        end
        L(a) = t ;
    end
    [newbestr,newbestant] = min(L) ;          %寻本次迭代最短路径及其相应蚂蚁
    if newbestr < bestr                       %到目前为止最优值的保存
        bestr = newbestr ;
        bestroad = visited(newbestant,:) ;
    end
    %离线更新信息素矩阵
    %挥发
    for a = 1 : m
        tau(visited(a,n),visited(a,1)) = tau(visited(a,n),visited(a,1)) * (1 - rho) ;
        for b = 1 : (n - 1)
            tau(visited(a,b),visited(a,(b + 1))) = tau(visited(a,b),visited(a,(b + 1))) * (1 - rho) ;
        end
    end
    %加强
    tau(visited(newbestant,n),visited(newbestant,1)) = tau(visited(newbestant,n),visited(newbestant,1)) + Q / L(newbestant) ;
    for b = 1 : (n - 1)
        tau(visited(newbestant,b),visited(newbestant,(b + 1))) = tau(visited(newbestant,b),visited(newbestant,(b + 1))) + Q / L(newbestant) ;
    end
end
bestr
bestx = zeros(1,n) ;
besty = zeros(1,n) ;
for i = 1 : n
    bestx(i) = x(bestroad(i)) ;
    besty(i) = y(bestroad(i)) ;
end
bestx = [bestx,bestx(1)] ;
besty = [besty,besty(1)] ;
plot(bestx,besty,'-') ;