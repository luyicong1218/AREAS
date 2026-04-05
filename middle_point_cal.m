function [middle_point_new] = middle_point_cal(max_radius,middle_point,middle_unit,member_new,point,vector)
%球节点最大半径，中点坐标，模型三维信息
middle_point_new = zeros(size(middle_point,1),3); %用于存储精确的节点坐标
member_all = []; %用于存储属于某个单元的所有节点，以便于减少对比的数量
member_ball_all = 1:1:size(point,1); %产生1到点云数量的数组
all_color = ['r','y','g','b','c','m','k']; %存储颜色字符用于不同杆件区别作图
for i = 1:length(member_new)
    member_all = [member_all,member_new{i}];
end
member_ball_all = setdiff(member_ball_all,member_all); %得到所有不属于杆件单元的节点集合
for i = 1:length(middle_unit)
    if size(middle_unit{i},2) == 0
        continue
    end
    if size(middle_unit{i},2) == 1
        middle_point_new(i,:) = middle_point(i,:);
        continue
    end
    %声明member_ball{i}用于存储属于第i个节点的所有点云
    member_ball{i} = [];
    for j = 1:size(member_ball_all,2)
        if norm(point(member_ball_all(j),:) - middle_point(i,:)) < 1.2 * max_radius && dot(point(member_ball_all(j),:) - middle_point(i,:),[vector.vertex.nx(member_ball_all(j)),vector.vertex.ny(member_ball_all(j)),vector.vertex.nz(member_ball_all(j))]) / norm(point(member_ball_all(j),:) - middle_point(i,:)) / norm([vector.vertex.nx(member_ball_all(j)),vector.vertex.ny(member_ball_all(j)),vector.vertex.nz(member_ball_all(j))]) > 0.85
            member_ball{i} = [member_ball{i}
                              point(member_ball_all(j),:)];
        else
            continue
        end
    end
    %进行球节点拟合
    %%球拟合，求出球心位置，及球的直径
    if size(member_ball{i},1) < 50
        middle_point_new(i,:) = middle_point(i,:);
        continue
    end
    data=unique(member_ball{i},'rows');
%     [middle_point_new(i,:),~] = sphereFit(data);
   f=@(p)(data(:,1)-p(1)).^2+(data(:,2)-p(2)).^2+(data(:,3)-p(3)).^2-p(4)^2;
   p=lsqnonlin(f,[middle_point(i,1) middle_point(i,2) middle_point(i,3) max_radius]');%拟合的参数
%    p=ransac(data,max_radius,middle_point(i,:));
   middle_point_new(i,:) = [p(1),p(2),p(3)];
end

fprintf('ball image outputing......\n');
figure('name','structure joint')
for i = 1:length(member_ball)
    if size(member_ball{i},2) == 0
        continue;
    end
    if rem(i,7) == 0
        color = all_color(7);
    else
        color = all_color(rem(i,7));
    end
    scatter3(member_ball{i}(:,1),member_ball{i}(:,2),member_ball{i}(:,3),'.',color);
    axis equal
    hold on
end

        