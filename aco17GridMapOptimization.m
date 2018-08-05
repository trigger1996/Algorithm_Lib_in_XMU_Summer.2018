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

G = ~G;

L = size(G, 1);

% ��ά����Ч������
%Weight_Mat = [ 10, 5,   2;
%			   5,  1,   0.5;
%			   2,  0.5, 0.1; ];

%G_conv2 = conv2(G, Weight_Mat);

startGrid = [ 0, 0 ];
endGrid   = [ L, L ];


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

% �����ڽӾ���
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
		   