%�����㷨test 
%�ò�����һ��Բ�ϵ�ʮ���������������㷨
  
clc
clear
  
%����
alpha = 1 ;                               %��Ϣ��ָ��
beta = 5  ;                                %����ָ��
rho = 0.5 ;                                %�ӷ�ϵ��
n = 16 ;                                   %���и���
k = 20 ;                                    %��������
m = n - 1 ;                                %����ֻ��������ȡ�ȳ�����Ŀ��һ������ֻ�� n - 1
Q = 100 ;
bestr = inf ;
%����һ��Բ�ϵ�ʮ����
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
%�������
d = zeros(n,n) ;
for i = 1 : n
    for j = 1 : n
        d(i,j) = sqrt( ( x(i) - x(j) ) ^ 2 + ( y(i) - y(j) ) ^ 2) ;
    end
end
temp = min(d) ;
dmin = temp(1) ;
tau = ones(n,n) ;
%tau = tau ./ (n * dmin) ;                   %��ʼ��tau��Ϣ�ؾ���
  
%��ʼ����
for i = 1 : k  
    %��ʼ��
    visited = zeros(m,n) ;                  %��visited ���������������߹��ĳ��� m��n ����δ����ĳ���Ϊ0
    visited(:,1) = (randperm(n,m))';        %��mֻ�����������n������ ������һ��1��n����������е�һ�����ݵĸ���
    for b = 2 : n                           %�������϶��ߵ���b������ʱ
        current = visited(:,(b-1)) ;         %���������������ڳ��� m��1
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
        for a = 1 : m                       %������aֻ����
            p = zeros(1,l) ;
            for j = 1 : l                   %������ʽ��ѡ����һ������
                p(j) = ( ( tau( current(a,1) , allow(a,j) ) ) ^ alpha ) * ( ( 1 / d( current(a,1) , allow(a,j) ) ) ^ beta ) ;
            end
            p = p ./ sum(p) ;               %�������̶ĵķ�ʽ
            p = cumsum(p) ;
            pick = rand ;
            for c = 1 : l
                if pick < p(c)
                    visited(a,b) = allow(a,c) ;          %�ҵ�����Ҫ��ĳ��� �� ��������a��·����
                    break ;
                end
            end
        end
    end
    %����ÿֻ�������ߵ�·���ܳ�
    L = zeros(1,m) ;
    for a = 1 : m
        t = d(visited(a,n),visited(a,1)) ;
        for b = 1 : (n - 1)
            t = t + d(visited(a,b),visited(a,(b + 1)));
        end
        L(a) = t ;
    end
    [newbestr,newbestant] = min(L) ;          %Ѱ���ε������·��������Ӧ����
    if newbestr < bestr                       %��ĿǰΪֹ����ֵ�ı���
        bestr = newbestr ;
        bestroad = visited(newbestant,:) ;
    end
    %���߸�����Ϣ�ؾ���
    %�ӷ�
    for a = 1 : m
        tau(visited(a,n),visited(a,1)) = tau(visited(a,n),visited(a,1)) * (1 - rho) ;
        for b = 1 : (n - 1)
            tau(visited(a,b),visited(a,(b + 1))) = tau(visited(a,b),visited(a,(b + 1))) * (1 - rho) ;
        end
    end
    %��ǿ
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