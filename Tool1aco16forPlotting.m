function Tool1aco16forPlotting(antTraj, AO_Length, G, a)

figure;

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
bestTraj = antTraj(1, :);		% ��ȡ��õ���ֻ���ϵ�·����˳������һ��0ȥ��
bestTraj = bestTraj(1 : size(bestTraj, 2) - 1);
bestTraj(:, find(sum(bestTraj, 1) == 0)) = [];						% ����ȫ0��

x = zeros(1, size(bestTraj, 2));
y = zeros(1, size(bestTraj, 2));
for i = 1 : size(bestTraj, 2)
	x(i)=a * (mod(bestTraj(i), MM) -0.5);
	y(i)=a * (MM + 0.5 - ceil(bestTraj(i) / MM));
end
plot(x, y);
end