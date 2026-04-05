function[id2] = endpoint_cal(member,num,vec,point)
%变量声明与初始化
ex_point = [];
%变量声明与初始化完成

%进行坐标变换，求解点云的分布长度，返回标记值
for i = 1:num    
    %对杆件的点云进行坐标变换，变换为与z轴的正方向一致
    id2{i} = [];
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
    if max(ex_point(:,3)) - min(ex_point(:,3)) < 100 %原来是50
        id2{i} = [id2{i},1,1];
    else
        id2{i} = [id2{i},1];
    end
    ex_point = [];
end
%求解点云分布长度，返回标记值完成