clc
clear 

%读取.ply型数据，根据位置与拓扑关系确定点云集团子集
%mesh存储点云的拓扑关系，point存储点云的几何位置坐标，vector存储顶点的法向量
fprintf('data reading......\n');
%%[mesh,point,vector] = plyread('E:\20190913 first paper\model\sanweisaomiao\wangjia2-1small.ply','tri');
[mesh,point,vector] = plyread('wangjia1-1.ply','tri');

while 1
    unit = input('please choose the unit of the data \n1. m\n2. dm\n3. cm\n4. mm\nplease input below number of the unit\nnumber:');
    if unit == 1 || unit == 2 || unit == 3 || unit == 4
        switch unit
            case 1
                point = 1000 * point;
            case 2
                point = 100 * point;
            case 3
                point = 10 * point;
            case 4
                point = 1 * point;
        end
        break;
    else
        fprintf('input error,please input again!\n');
    end
end
max_radius = input('please input the max radius of the ball:');
max_bending = input('please input the max ratio of deflection to span of the bar:');
fprintf('data read finished.\n');
%.ply数据读取完毕

%变量声明与初始化
points = [transpose([1:size(point,1)]),point];%points为将点云位置坐标数组添加节点编号索引后的数组
node = []; %初始结点集合
num = 1;%杆件编号初始值
while_num = 0;%指示while循环的进行次数
mesh_row = [];%mesh中与索引对应的所有行号（之后要对其中所有行删除）
points_row = [];%points中与索引对应的所有行号（之后要对其中所有行删除）
node_num = 0;%计算在单个索引中属于节点的点的个数
all_color = ['r','y','g','b','c','m','k']; %存储颜色字符用于不同杆件区别作图
%对数据进行筛除，删去法向量信息无效的点，即x、y、z方向法向量分量为0的点以及该点与其他点之间的联系
angle_sum = [];%对所有判断是否为单元的角值汇总，以便于设置合适的阈值
angle_node = [];%对所有判断是否为单元的角值汇总，以便于设置合适的阈值
vec = [];%存储所有杆件的主向量信息
mesh_cpy = mesh;%存储mesh的初始集合
%变量声明与初始化完毕

%对法向量为0的无用点进行删除
[row,~] = find(vector.vertex.nx == 0);
for i = 1:size(row,1)
    [row1,~] = find(mesh == row(i));
    mesh_row = unique([mesh_row,transpose(row1)]);
end
points(row,:) = [];
mesh(mesh_row,:) = [];
mesh_row = [];
meshes = mesh; %将mesh数据赋值给meshes，保留所用有用的mesh数据
%对法向量为0的无用点删除完成

%进行单元的识别与划分，从第一个点开始寻找关联点
fprintf('member dividing (step 1)......\n');
while (size(points,1) ~= 0)
    index = [points(1,1)];%添加索引的初始值
    member{num} = [];
    angle_node = [];
    while size(index,2) ~= 0
        while_num = while_num + 1;  %while循环次数递增        
        points_row = index;
        %得到该点的下一层所有结点的索引值
        for i = 1:size(index,2)
            [row,col] = find(mesh == index(i));
            mesh_row = unique([mesh_row,transpose(row)]);
            for j =1:size(row,1) - isempty(row)
                switch col(j)
                    case 1
                        index = [index,mesh(row(j),2),mesh(row(j),3)];
                    case 2
                        index = [index,mesh(row(j),1),mesh(row(j),3)];
                    case 3
                        index = [index,mesh(row(j),1),mesh(row(j),2)];
                end
            end
        end
        index(1:size(points_row,2)) = []; %删除初始赋值 
        index = unique(index); %得到与该节点相关的无重复的其他结点编号
        %在第一次进入循环后需要计算主向量，判断该点是否为节点，如果是则跳出循环，如果不是则归入杆件点云集合
        if while_num == 1
            if size(index,2) == 0
                node = [node,points(1,1)];
                points(1,:) = [];
                break;
            end
            [angle_sum{num},p_vector,judge] = angle_judgement(points(1,1),index,vector);
            if judge == 0
                member{num} = [member{num},points(1,1),index]; %将初始点置入第num个杆件点云集合
                vec = [vec,transpose(p_vector)]; %将该杆件的主向量存储在vec数组中，以便于后续结点的判断
                mesh(mesh_row,:) = [];
                points(1,:) = [];
                mesh_row =[];
                continue;
            else
                node = [node,points(1,1)];
                mesh(mesh_row,:) = [];
                mesh_row = [];
                points(1,:) = [];
                break;
            end
        end         
        %当while循环进行第二次时开始执行下面的程序
        loop_num = size(index,2); %计算循环次数
        for i = 1:loop_num
             if abs(dot(p_vector,[vector.vertex.nx(index(i - node_num)),vector.vertex.ny(index(i - node_num)),vector.vertex.nz(index(i - node_num))]) / (norm(p_vector) * norm([vector.vertex.nx(index(i - node_num)),vector.vertex.ny(index(i - node_num)),vector.vertex.nz(index(i - node_num))]))) < 0.2 && member_judgement(p_vector,index(i - node_num),meshes,vector) == 0
                   member{num} = [member{num},index(i - node_num)];
             else
                   node = [node,index(i - node_num)];
                   [n_row,n_col] = find(mesh == index(i -node_num));
                   mesh_row = unique([mesh_row,transpose(n_row)]);
                   index(i - node_num) = [];
                   node_num = node_num + 1;                  
             end
        end
        mesh(mesh_row,:) = [];
        for k =1:size(points_row,2)    
            points(find(points(:,1) == points_row(k)),:) = [];
        end
        mesh_row =[];
        node_num = 0;   
    end
    %判断是否有数据存入定义的点云集合中，如果没有，仍然将后续数据存在该集合，如果有，自动生成下一个细胞数组 
    if isempty(member{num}) ~= 1
        num = num + 1; %num表示第num个杆件的点云集合，用于细胞数组的存储
    end
    while_num = 0;
end
fprintf('member divide (step 1) finished.\n');
%单元的粗略识别与划分完成

%进行单元集合的精细筛除，除去不符合要求的单元集合
 %针对所有杆件中点云数量不大于1/20平均值的元胞数组进行删除，认为其不构成一个单元
fprintf('member dividing (step 2)......\n');
id1 = cellfun('length',member) < max(cellfun('length',member)) / 50;
member(id1) = [];
if size(id1,2) > size(vec,2)
    vec(:,setdiff(find(id1 == 1),[size(id1,2)])) = [];
else
    vec(:,find(id1 == 1)) = [];
end
num = length(member);%更新杆件数量
 %其次利用单元的长度不大于规定值判据将单元信息进行更新，删除不符合要求的元胞数组
id2 = cellfun('length',endpoint_cal(member,num,vec,point)) == 2;
member(id2) = [];
vec(:,find(id2 == 1)) = [];
num = length(member);%再次更新杆件数量
fprintf('member divide (step 2) finished.\n');
%单元集合精细筛除完毕    
%进行节点的判断，判断的过程中对两端点相同的杆件信息进行合并，认为其属于同一个单元
fprintf('node location calculating......\n');
[middle_unit,middle_point,member_index] = node_location(member,num,vec,point,vector); %返回值为节点的集合值
fprintf('node location calculated.\n');
%节点判断完成

%进行杆件变形的识别与提取
fprintf('member movement calculating......\n');
member_new = member_div(max_radius,max_bending,middle_unit,middle_point,member_index,point,num,vec,vector,mesh_cpy);
fprintf('member movement calculated.\n');
%杆件变形的识别与提取完成

%将所有的节点作图显示
fprintf('image outputing......\n');
figure('name','structure model')
% scatter3(middle_point(:,1),middle_point(:,2),middle_point(:,3),'*','r')
% axis equal
% hold on
for i = 1:num
    if size(member_new{i},2) == 0
        continue;
    end
    if rem(i,7) == 0
        color = all_color(7);
    else
        color = all_color(rem(i,7));
    end
    scatter3(point(transpose(member_new{i}),1),point(transpose(member_new{i}),2),point(transpose(member_new{i}),3),'.',color)
    axis equal
    text(point(member_new{i}(1),1),point(member_new{i}(1),2),point(member_new{i}(1),3),num2str(i)) %将杆件编号在图中显示
    hold on
end
fprintf('image output finished.\n');
%单元点云画图显示完成     

    




