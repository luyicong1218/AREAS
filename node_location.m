function[middle_unit,middle_point,member_index] = node_location(member,num,vec,point,vector)
%变量声明与初始化
ex_point = [];%该变量存储变换之后的点云x、y、z坐标
index = []; %index存储某个杆件点云在变换之后z坐标的最大最小值所在的位置索引，按照单元顺序进行存储
line = [];%存储index值的行向量
col = []; %记录要删除的行号索引值
member_index = zeros(num,2); %记录单元所对应的节点编号索引值
r = [];%用于存储所得的杆件直径半径数据，以便于对节点的位置坐标进行修正
x = [];
y = [];
%变量声明与初始化完成

%求取每个杆件单元的两端点索引值
%对所有杆件单元的点云做如下操作：将杆件的主向量方向百年换至与坐标轴z方向一致，之后寻找z坐标的最大最小值以确定杆件的两端点
for i = 1:num    
    %对杆件的点云进行坐标变换，变换为与z轴的正方向一致
    node_vector = transpose(vec(:,i));%将单元的主向量值赋给vector变量进行计算
    %求解变换矩阵
    Mat1 = [1,0,0,0
            0,node_vector(3) / sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2),node_vector(2) / sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2),0
            0,-node_vector(2) / sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2),node_vector(3) / sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2),0
            0,0,0,1];
    Mat2 = [sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2) / norm(node_vector),0,node_vector(1) / norm(node_vector),0
            0,1,0,0
            -node_vector(1) / norm(node_vector),0,sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2) / norm(node_vector),0
            0,0,0,1];
    Mat = Mat1 * Mat2;
    for j = 1:size(member{i},2)
        ex_point = [ex_point
                    [point(member{i}(j),:),1] * Mat];
    end
    [row1,~] = find(ex_point(:,3) == max(ex_point(:,3)));
    [row2,~] = find(ex_point(:,3) == min(ex_point(:,3)));
    index(i,:) = [member{i}(row2(1)),member{i}(row1(1))];
    ex_point = [];
end
%求取每个杆件单元两端点索引值完毕

%求取节点点云集合以及杆件所对应的的节点编号索引
num = 0 ;%将num值赋0，用于统计节点的个数
for i = 1:size(index,1)
    line = [line,index(i,:)];%将index的值赋给line向量
end
while numel(line) ~= 0
    num = num + 1;
    node{num} = [];
    for j = 1:size(line,2)
        if sqrt((point(line(1),1) - point(line(j),1)) ^ 2 + (point(line(1),2) - point(line(j),2)) ^ 2 + (point(line(1),3) - point(line(j),3)) ^ 2) < 100 %设定点之间的阈值为800，该值与球节点的直径大小有关
            node{num} = [node{num},line(j)];
            col = [col ,j];
            [row,~] = find(index == line(j));
            col1 = find(member_index(row(1),:) == 0);
            member_index(row(1),col1(1)) = num;
        end
    end
    node{num} = unique(node{num});
    line(transpose(col)) = [];
    col = [];
end
for i = 1:size(member_index,1)
    member_index(i,:) = sort(member_index(i,:));%对编号从小到大排列
end
%节点点云集合以及杆件所对应的的节点编号索引求取完毕

%求取杆件单元的半径大小
for i = 1:length(member)
    %对杆件的点云进行坐标变换，变换为与z轴的正方向一致
    node_vector = transpose(vec(:,i));%将单元的主向量值赋给vector变量进行计算
    %求解变换矩阵
    Mat1 = [1,0,0,0
            0,node_vector(3) / sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2),node_vector(2) / sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2),0
            0,-node_vector(2) / sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2),node_vector(3) / sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2),0
            0,0,0,1];
    Mat2 = [sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2) / norm(node_vector),0,node_vector(1) / norm(node_vector),0
            0,1,0,0
            -node_vector(1) / norm(node_vector),0,sqrt(node_vector(2) ^ 2 + node_vector(3) ^ 2) / norm(node_vector),0
            0,0,0,1];
    Mat = Mat1 * Mat2;
    for j = 1:size(member{i},2)
        ex_point = [ex_point
                    [point(member{i}(j),:),1] * Mat];
    end
    [row,~] = find(ex_point(:,3) > max(ex_point(:,3)) / 2 + min(ex_point(:,3)) / 2 - 50 & ex_point(:,3) < max(ex_point(:,3)) / 2 + min(ex_point(:,3)) / 2 + 50);
    x = [x,ex_point(row,1)];
    y = [y,ex_point(row,2)];
    circle_radius = CircleFitByTaubin([x,y]);  
    r = [r,circle_radius(3)];
    ex_point = [];
    x = [];
    y = [];
end
%杆件单元半径求取完毕
%对杆件单元半径进行修正
for i = 1:size(r,2)
    if r(i) > 100 * r(1)
        r(i) = r(i - 1);
    end
end
%对杆件单元半径修正完成

%同一节点的点云归一化处理
middle_point = zeros(num,3);%用于存储归一化后节点的三维坐标
%middle_unit{i}用于存储第i个节点所连接的所有单元的编号
for i = 1:num
    if size(node{i},2) == 1
        [row,~] = find(index == node{i}(1));
        middle_point(i,:) = point(node{i}(1),:) - r(1,row(1)) * [vector.vertex.nx(node{i}(1)),vector.vertex.ny(node{i}(1)),vector.vertex.nz(node{i}(1))] / norm([vector.vertex.nx(node{i}(1)),vector.vertex.ny(node{i}(1)),vector.vertex.nz(node{i}(1))]);
        middle_unit{i} = [row(1)];
        continue;
    end
    middle_unit{i} = [];
    for j = 1:size(node{i},2) - 1
        [row,~] = find(index == node{i}(j));
        vec1 = vec(:,row);
        middle_unit{i} = [middle_unit{i},row(1)];
        [row,~] = find(index == node{i}(j + 1));
        vec2 = vec(:,row);
        middle_unit{i} = [middle_unit{i},row(1)];
        [row,~] = find(index == node{i}(j));
        [row1,~] = find(index == node{i}(j + 1));
        if abs(dot(vec1,vec2) / (norm(vec1) * norm(vec2))) > 0.9
            middle_point(i,:) = middle_point(i,:) + (point(node{i}(j),:) - r(row(1)) * [vector.vertex.nx(node{i}(j)),vector.vertex.ny(node{i}(j)),vector.vertex.nz(node{i}(j))] / norm([vector.vertex.nx(node{i}(j)),vector.vertex.ny(node{i}(j)),vector.vertex.nz(node{i}(j))]) + point(node{i}(j + 1),:) - r(row1(1)) * [vector.vertex.nx(node{i}(j + 1)),vector.vertex.ny(node{i}(j + 1)),vector.vertex.nz(node{i}(j + 1))] / norm([vector.vertex.nx(node{i}(j + 1)),vector.vertex.ny(node{i}(j + 1)),vector.vertex.nz(node{i}(j + 1))])) / 2;
        else
            [~,mid] = line_middle_point(point(node{i}(j),:) - r(1,row(1)) * [vector.vertex.nx(node{i}(j)),vector.vertex.ny(node{i}(j)),vector.vertex.nz(node{i}(j))] / norm([vector.vertex.nx(node{i}(j)),vector.vertex.ny(node{i}(j)),vector.vertex.nz(node{i}(j))]),vec1,point(node{i}(j + 1),:) - r(row1(1)) * [vector.vertex.nx(node{i}(j + 1)),vector.vertex.ny(node{i}(j + 1)),vector.vertex.nz(node{i}(j + 1))] / norm([vector.vertex.nx(node{i}(j + 1)),vector.vertex.ny(node{i}(j + 1)),vector.vertex.nz(node{i}(j + 1))]),vec2);
            middle_point(i,:) = middle_point(i,:) + mid;
        end
    end
    middle_unit{i} = unique(middle_unit{i});
    middle_point(i,:) = middle_point(i,:) / (size(node{i},2) - 1);
end
%点云归一化处理完毕

