clear;
clc;
clear all;

f = @(x) x .* sin(10 * pi * x) + 2.;% f(x) = sin(10 * pi * x) + 2, x取值范围[-1, 2]
%f = @(x) sin(x) + x .* cos(x);		% 函数表达式
fplot(f, [-1, 2]);					% 画出函数图像

N = 50;								% 种群上限
L = 10;								% 基因长度
iter = 100;							% 迭代次数
dna = randi([0, 1], [N, L]);		% 基因，rand()参数1是每个数的取值范围，参数2是矩阵的行列值

% 变异
% 变异最简单，就是个体中的某一位的值发生了变化
x3 = dna;
for i = 1: N                        % 变异操作
    if rand < 1						% 本来变异是个几率，而且很低，但是这里为了做实验，强制变异
        x3(i,randi(L)) = randi([0, 1]);
    end
end