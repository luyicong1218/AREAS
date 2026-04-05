function [results] = member_judgement(vec,point,mesh,vector)
%变量声明与初始化
temp = 0;
index = []; %存储与该点有关系的其他所有点编号
%变量声明与初始化完成

%主向量求取
[row,col] = find(mesh == point);
for j =1:size(row,1)
    switch col(j)
    case 1
        index = [index,mesh(row(j),2),mesh(row(j),3)];
    case 2
        index = [index,mesh(row(j),1),mesh(row(j),3)];
    case 3
        index = [index,mesh(row(j),1),mesh(row(j),2)];
    end
end
index = unique(index); %得到与该节点相关的无重复的其他结点编号
results = 0; %存储判断参数，为0表示满足法向量共面，否则表示法向量不共面，该节点应该为结构中的结点
if size(index,2) < 2
    results = 1;
else
    p_vector = zeros(size(index,2),3); %存储每个主向量的x、y、z坐标
    p_vector(1,:) = cross([vector.vertex.nx(index(1)),vector.vertex.ny(index(1)),vector.vertex.nz(index(1))],[vector.vertex.nx(index(size(index,2))),vector.vertex.ny(index(size(index,2))),vector.vertex.nz(index(size(index,2)))]);
    for i = 1:size(index,2) - 1
        p_vector(i + 1,:) = cross([vector.vertex.nx(index(i)),vector.vertex.ny(index(i)),vector.vertex.nz(index(i))],[vector.vertex.nx(index(i + 1)),vector.vertex.ny(index(i + 1)),vector.vertex.nz(index(i + 1))]);
    end
    %主向量求取完成

    %判断主向量是否共线
    for i = 1:size(index,2)
        if abs(dot(vec,p_vector(i,:))/(norm(vec) * norm(p_vector(i,:)))) < 0.9
            temp = temp + 1;
        end
    end
    if temp >= ceil(size(index,2) / 2)
        results = 1;
    end
end
%判断主向量是否共线完成
