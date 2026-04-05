function [angle_sum,vec,result] = angle_judgement(point,index,vector)
%变量声明与初始化
angle_sum = [];
vec = zeros(1,3); %存储主向量坐标
result = 0; %存储判断参数，为0表示满足法向量共面，否则表示法向量不共面，该节点应该为结构中的结点
%变量声明与初始化完成

%求解主向量
p_vector = zeros(size(index,2),3); %存储每个主向量的x、y、z坐标
p_vector(1,:) = cross([vector.vertex.nx(point),vector.vertex.ny(point),vector.vertex.nz(point)],[vector.vertex.nx(index(1)),vector.vertex.ny(index(1)),vector.vertex.nz(index(1))]);
for i = 1:size(index,2) - 1
    p_vector(i + 1,:) = cross([vector.vertex.nx(index(i)),vector.vertex.ny(index(i)),vector.vertex.nz(index(i))],[vector.vertex.nx(index(i + 1)),vector.vertex.ny(index(i + 1)),vector.vertex.nz(index(i + 1))]);
end
%主向量求解完成

%判断所有的主向量是否满足共线条件
for i = 1:size(index,2) - 1
    angle_sum = [angle_sum,abs(dot(p_vector(1,:),p_vector((i + 1),:))/(norm(p_vector(1,:)) * norm(p_vector((i + 1),:))))];
    if abs(dot(p_vector(1,:),p_vector((i + 1),:))/(norm(p_vector(1,:)) * norm(p_vector((i + 1),:)))) < 0.95   %多视角三维重建是0.9
        result = 1;
    end
end
%判断完成

%将所有主向量求和平均值作为返回的主向量值
for i = 1:size(index,2)
    vec = vec + sign(p_vector(1,1) * p_vector(i,1)) * p_vector(i,:) / (norm(p_vector(i,:)) * size(index,2));
end
%求解主向量最终值完成
        
