function[member_new] = member_div(max_radius,max_bending,middle_unit,middle_point,member_index,point,num,vec,vector,mesh_cpy)
%变量声明与初始化
del = [];%用于存储单元间主向量的数量积
del_unit = [];%记录需要删除的单元索引
ex_point = [];%存储坐标变换后的三维坐标
member_point = [];%存储20个基准点的坐标
% logo_point = [];%用于记录各个标记点修正后的坐标
distance = 0;%记录变形值
judgement = 0;%变形的判定标记
x = [];
y = [];
distance2 = [];
%变量声明与初始化完成

%进行点云的后续操作，如由于遮挡引起的杆件断开，实则属于同一单元
for i = 1:size(member_index,1)
%若一个单元对应的两个节点索引相同，则该单元无效，将其单元对应的节点索引置零，同时，将节点对应的单元索引中除去该单元
    if member_index(i,1) == member_index(i,2) 
        middle_unit{member_index(i,1)} = setdiff(middle_unit{member_index(i,1)},i);
        member_index(i,:) = zeros(1,2);
    end
end

for i = 1:length(middle_unit)
    %此处对应的是悬臂杆件，即杆件的一端不与任何其他杆件通过球节点相连，这种情况下不需要进行单元的合并等操作，直接跳过
    if size(middle_unit{i},2) < 2
        continue;
    end
    %对于非悬臂杆件，判断需要删除的杆件，若两个单元共线，则该单元无效，将其单元对应的节点索引置零，同时，将节点对应的单元索引中除去该单元
    for j = 1:size(middle_unit{i},2) - 1
        for k = j + 1:size(middle_unit{i},2)
            nod1 = setdiff(member_index(middle_unit{i}(j),:),i);
            vec1 = middle_point(nod1(1),:) - middle_point(i,:);
            nod2 = setdiff(member_index(middle_unit{i}(k),:),i);
            vec2 = middle_point(nod2(1),:) - middle_point(i,:);
            del = [del
                   middle_unit{i}(j),middle_unit{i}(k),nod1(1),nod2(1),vec1,vec2,dot(vec1,vec2) / norm(vec1) / norm(vec2)];
        end
    end
    [row,~] = find(del(:,11) > 0.8);
    if ~isempty(row)
        for j = 1:size(row,1)
            if ismember(del(row(j),1),del_unit) == 1 || ismember(del(row(j),2),del_unit) == 1
                continue;
            end
            switch find([norm(del(row(j),5)),norm(del(row(j),6))] == max(norm(del(row(j),5)),norm(del(row(j),6))))
                case 1
                    middle_unit{del(row(j),4)} = setdiff(middle_unit{del(row(j),4)},del(row(j),2));
                    member_index(del(row(j),2),:) = zeros(1,2);
                    vec(:,del(row(j),2)) = zeros(3,1);
                    del_unit = [del_unit,del(row(j),2)];
                case 2
                    middle_unit{del(row(j),3)} = setdiff(middle_unit{del(row(j),3)},del(row(j),1));
                    member_index(del(row(j),1),:) = zeros(1,2);
                    vec(:,del(row(j),1)) = zeros(3,1);
                    del_unit = [del_unit,del(row(j),1)];
            end
        end
    end
    [row,~] = find(del(:,11) < -0.8);
    if ~isempty(row) && size(middle_unit{i},2) - size(del_unit,2) == 2
        for j = 1:size(row,1)
            if ismember(del(row(j),1),del_unit) == 1 || ismember(del(row(j),2),del_unit) == 1
                continue;
            end
            member_index(del(row(j),2),:) = sort([del(row(j),3),del(row(j),4)]);
            member_index(del(row(j),1),:) = zeros(1,2);
            vec(:,del(row(j),1)) = zeros(3,1);
            middle_unit{del(row(j),3)} = [setdiff(middle_unit{del(row(j),3)},del(row(j),1)),del(row(j),2)];
            middle_unit{i} = [];
        end
    end
    for j = 1:size(del_unit,2)
        middle_unit{i} = setdiff(middle_unit{i},del_unit(j));
    end
    del_unit = [];
    del = [];
end
%点云后续操作完成

%进行点云的重新划分，此次划分利用位置关系进行判定，不依赖于拓扑关系，因此点云密度很高
for i = 1:num
    member_new{i} = [];
    if member_index(i,:) == zeros(1,2)
        continue;
    end
    judgement = [abs(dot(transpose(vec(:,i)),[1,0,0]) / norm(vec(:,i))),abs(dot(transpose(vec(:,i)),[0,1,0]) / norm(vec(:,i))),abs(dot(transpose(vec(:,i)),[0,0,1]) / norm(vec(:,i)))];
    col = find(judgement == max(judgement));
    %对中点坐标进行修正，除去球节点部分的变形
    node1 = middle_point(member_index(i,1),:) + (middle_point(member_index(i,2),:) - middle_point(member_index(i,1),:)) / norm([middle_point(member_index(i,2),:) - middle_point(member_index(i,1),:)]) * max_radius;
    node2 = middle_point(member_index(i,2),:) + (middle_point(member_index(i,1),:) - middle_point(member_index(i,2),:)) / norm([middle_point(member_index(i,2),:) - middle_point(member_index(i,1),:)]) * max_radius;
    row = find(point(:,col(1)) > min(node1(col(1)),node2(col(1))) & point(:,col(1)) < max(node1(col(1)),node2(col(1))));
    for j = 1:size(row,1)
        temp = middle_point(member_index(i,2),:) - middle_point(member_index(i,1),:);
        if norm(cross(point(row(j),:) - middle_point(member_index(i,1),:),middle_point(member_index(i,1),:) - middle_point(member_index(i,2),:))) / norm(middle_point(member_index(i,1),:) - middle_point(member_index(i,2),:)) < 200 * sin (dot([point(row(j),:) - middle_point(member_index(i,1),:)],temp / norm(temp)) / norm(temp) * pi) && member_judgement(transpose(vec(:,i)),row(j),mesh_cpy,vector) == 0
            member_new{i} = [member_new{i},row(j)];
        end
    end
end

%进行一次节点的精确定位
middle_point_new = middle_point_cal(max_radius,middle_point,middle_unit,member_new,point,vector);
%节点精确定位完成
fprintf('possible bending members as follows.\n');
for i = 1:num
    if isempty(member_new{i}) == 1
        continue;
    end
    vector1 = middle_point_new(member_index(i,2),:) - middle_point_new(member_index(i,1),:);    %利用所得的结点求该杆件的主向量
    %利用主向量进行坐标变换，使得主向量方向与z轴方向相同
    Mat1 = [1,0,0,0
            0,vector1(3) / sqrt(vector1(2) ^ 2 + vector1(3) ^ 2),vector1(2) / sqrt(vector1(2) ^ 2 + vector1(3) ^ 2),0
            0,-vector1(2) / sqrt(vector1(2) ^ 2 + vector1(3) ^ 2),vector1(3) / sqrt(vector1(2) ^ 2 + vector1(3) ^ 2),0
            0,0,0,1];
    Mat2 = [sqrt(vector1(2) ^ 2 + vector1(3) ^ 2) / norm(vector1),0,vector1(1) / norm(vector1),0
            0,1,0,0
            -vector1(1) / norm(vector1),0,sqrt(vector1(2) ^ 2 + vector1(3) ^ 2) / norm(vector1),0
            0,0,0,1];
    Mat = Mat1 * Mat2;
    for j = 1:size(member_new{i},2)
        ex_point = [ex_point
                    [point(member_new{i}(j),:),1] * Mat];
    end
    %ex_point数组存储所有属于该杆件单元的变换后三维坐标
    n_point =[[middle_point_new(member_index(i,1),:),1] * Mat         %middle_point_new
               [middle_point_new(member_index(i,2),:),1] * Mat];    
    %求取控制点，控制点的间距为50mm
    %求取控制点之前先对杆件的直径进行统一,防止出现过大的直径，拟合失效，变形计算出现偏差
    for j = 1:floor((max(ex_point(:,3)) - min(ex_point(:,3))) / 30)
        [row,~] = find(ex_point(:,3) >= min(ex_point(:,3)) + (j - 1) * 30 & ex_point(:,3) <= min(ex_point(:,3)) + j * 30);
        if size(row,1) < 10
            x=[];
            y=[];
            continue
        else
            x = [x,ex_point(row,1)];
            y = [y,ex_point(row,2)];
            circle_radius = CircleFitByTaubin([x,y]);
            f=@(p)(transpose(x) - p(1)).^2 + (transpose(y) - p(2)).^2 - 12.5 ^ 2;
%             p=ransac(transpose([x,y]),circle_radius);
            p=lsqnonlin(f,[circle_radius(1) circle_radius(2)]');%拟合的参数         
            x = [];
            y = [];
            member_point = [member_point
                            p(1),p(2),min(ex_point(:,3)) + (j - 1 / 2) * 30];
%             circle_radius = CircleFitByTaubin([x,y]);         
%             x = [];
%             y = [];
%             member_point = [member_point
%                             circle_radius(1),circle_radius(2),min(ex_point(:,3)) + (j - 1 / 2) * 15,circle_radius(3)];

        end
    end
%     delete_row = [];
%     if isempty(member_point) == 1
%         fprintf('member %d has no information about movement\n',i')
%         continue
%     end
%     ave = mean(member_point(:,4));
%     sve = std(member_point(:,4));
%     for k = 1:size(member_point,1)
%         if abs(member_point(k,4) - ave) - 1 * sve > 0
%             delete_row = [delete_row
%                           k];
%         else
%             continue
%         end
%     end
%     member_point(delete_row,:) = [];
%     
    %控制点求取完成，存储在member_point中
    
    %对得到的控制点进行处理，判断杆件的变形状况，并提取杆件变形值
    for j = 2:size(member_point,1) - 2
%          if norm(cross(member_point(j,1:3) - n_point(1,1:3),n_point(2,1:3) - n_point(1,1:3))) / norm(n_point(2,1:3) - n_point(1,1:3)) > 1 / max_bending * norm(n_point(2,1:3) - n_point(1,1:3)) && norm(cross(member_point(j + 1,1:3) - n_point(1,1:3),n_point(2,1:3) - n_point(1,1:3))) / norm(n_point(2,1:3) - n_point(1,1:3)) > 1 / max_bending * norm(n_point(2,1:3) - n_point(1,1:3)) && norm(cross(member_point(j + 2,1:3) - n_point(1,1:3),n_point(2,1:3) - n_point(1,1:3))) / norm(n_point(2,1:3) - n_point(1,1:3)) > 1 / max_bending * norm(n_point(2,1:3) - n_point(1,1:3));
         if norm(cross(member_point(j,1:3) - member_point(2,1:3),member_point(size(member_point,1)-1,1:3) - member_point(2,1:3))) / norm(member_point(size(member_point,1)-1,1:3) - member_point(2,1:3)) > 1 / max_bending * norm(n_point(2,1:3) - n_point(1,1:3));
%              distance = max([distance,norm(cross(member_point(j,1:3) - n_point(1,1:3),n_point(2,1:3) - n_point(1,1:3))) / norm(n_point(2,1:3) - n_point(1,1:3)),norm(cross(member_point(j + 1,1:3) - n_point(1,1:3),n_point(2,1:3) - n_point(1,1:3))) / norm(n_point(2,1:3) - n_point(1,1:3)),norm(cross(member_point(j + 2,1:3) - n_point(1,1:3),n_point(2,1:3) - n_point(1,1:3))) / norm(n_point(2,1:3) - n_point(1,1:3))]);
             distance = max([distance,norm(cross(member_point(j,1:3) - member_point(2,1:3),member_point(size(member_point,1)-1,1:3) - member_point(2,1:3))) / norm(member_point(size(member_point,1)-1,1:3) - member_point(2,1:3))]);
             judgement = 1;
         end
    end
    while judgement == 1
        fprintf('member %d :\n',i);
        fprintf('the location of member %d is between [%f,%f,%f] and [%f,%f,%f]\n',i,middle_point(member_index(i,1),1),middle_point(member_index(i,1),2),middle_point(member_index(i,1),3),middle_point(member_index(i,2),1),middle_point(member_index(i,2),2),middle_point(member_index(i,2),3));
        fprintf('the biggest movement is  %f \n',distance);
        figure('name','member movement')
        scatter3(member_point(:,1),member_point(:,2),member_point(:,3),'*','r')
        axis equal
        hold on
        scatter3(ex_point(:,1),ex_point(:,2),ex_point(:,3),'.','b')
        axis equal
        hold on
        judgement = 0;
        distance2 = [distance2
                     n_point(1,3),0];
        for j = 1:size(member_point,1)
            distance2 = [distance2
                         member_point(j,3),norm(cross(member_point(j,1:3) - n_point(1,1:3),n_point(2,1:3) - n_point(1,1:3))) / norm(n_point(2,1:3) - n_point(1,1:3))];       
        end
        distance2 = [distance2
                     n_point(2,3),0];
        distance2
    end
    distance = 0;
    distance2 = [];
    ex_point = [];
    member_point = [];
end

        

