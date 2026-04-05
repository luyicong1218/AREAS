function [distance,p] = line_middle_point(node1,vec1,node2,vec2)
%变量声明与初始化
temp = zeros(1,3);  %存储辅助直线与平面的交点坐标
c_node1 = zeros(1,3);%存储第一个最近距离点的坐标
c_node2 = zeros(1,3);%存储第二个最近距离点的坐标
%变量声明与初始化完成

%求解坐标
N = cross(vec1,vec2);
temp(1) = (N(1) ^ 2 * node2(1) + N(1) * N(2) * node2(2) + N(1) * N(3) * node2(3) - N(1) * N(2) * node1(2) - N(1) * N(3) * node1(3) + N(2) ^ 2 * node1(1) + N(3) ^ 2 * node1(1)) / (N(1) ^ 2 + N(2) ^ 2 + N(3) ^ 2);
temp(2) = (N(1) * N(2) * node2(1) + N(2) ^ 2 * node2(2) + N(2) * N(3) * node2(3) - N(1) * N(2) * node1(1) - N(2) * N(3) * node1(3) + N(3) ^ 2 * node1(2) + N(1) ^ 2 * node1(2)) / (N(1) ^ 2 + N(2) ^ 2 + N(3) ^ 2);
temp(3) = (N(1) * N(3) * node2(1) + N(2) * N(3) * node2(2) + N(3) ^ 2 * node2(3) - N(2) * N(3) * node1(2) - N(1) * N(3) * node1(1) + N(1) ^ 2 * node1(3) + N(2) ^ 2 * node1(3)) / (N(1) ^ 2 + N(2) ^ 2 + N(3) ^ 2);
c_node1(1) = (vec1(1) * node2(2) - vec1(1) * temp(2) - vec1(1) * vec2(2) / vec2(1) * node2(1) + vec1(2) * temp(1)) / (vec1(2) - vec1(1) * vec2(2) / vec2(1));
c_node1(2) = vec2(2) / vec2(1) * (c_node1(1) - node2(1)) + node2(2);
c_node1(3) = vec2(3) / vec2(1) * (c_node1(1) - node2(1)) + node2(3);
c_node2(1) = (N(1) * node1(2) - N(1) * c_node1(2) - vec1(2) / vec1(1) * N(1) * node1(1) + N(2) * c_node1(1)) / (N(2) - N(1) *vec1(2) / vec1(1));
c_node2(2) = vec1(2) / vec1(1) * (c_node2(1) - node1(1)) + node1(2);
c_node2(3) = vec1(3) / vec1(1) * (c_node2(1) - node1(1)) + node1(3);
distance = norm(c_node1 - c_node2);
p = (c_node1 + c_node2) / 2;
%求解坐标完成