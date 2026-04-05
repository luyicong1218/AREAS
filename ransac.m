function p = ransac(member_ball,ball_radius,ball)
Data = transpose(member_ball);
%%  参数初始化
nSampLen = 4;               %设定模型所依据的点数
nDataLen = size(Data, 2);   %数据长度
nIter = 100;                 %最大循环次数
dThreshold = 1;             %残差阈值
nMaxInlyerCount=-1;         %点数下限
A=zeros([3 1]);
B=zeros([3 1]);
C=zeros([3 1]);
D=zeros([3 1]);
%%  主循环
for i = 1:nIter 
   SampleMask = zeros([1 nDataLen]);   
    while sum( SampleMask ) ~= nSampLen
        ind = ceil(nDataLen .* rand(1, nSampLen - sum(SampleMask))); %抽样，选取nSampLen个不同的点
        SampleMask(ind) = 1;
    end    
    Sample = find( SampleMask );                                    %找出非零元素的索引值，即建立模型的点
    %%  建立模型，存储建模需要的圆的三个点
      %圆定义方程：到定点距离为常数
      A(:,1)=Data(:,Sample(1));    %圆上一点
      B(:,1)=Data(:,Sample(2));    %圆上一点     
      C(:,1)=Data(:,Sample(3));    %圆上一点
      D(:,1)=Data(:,Sample(4));
      mat=[transpose([A B C D]),transpose([1,1,1,1])];
      mat1=[transpose([-sum(A.^2),-sum(B.^2),-sum(C.^2),-sum(D.^2)]),mat(:,2:4)];
      mat2=[mat(:,1),transpose([-sum(A.^2),-sum(B.^2),-sum(C.^2),-sum(D.^2)]),mat(:,3:4)];
      mat3=[mat(:,1:2),transpose([-sum(A.^2),-sum(B.^2),-sum(C.^2),-sum(D.^2)]),mat(:,4)];
      mat4=[mat(:,1:3),transpose([-sum(A.^2),-sum(B.^2),-sum(C.^2),-sum(D.^2)])];

      temp_a=-det(mat1)/det(mat);
      temp_b=-det(mat2)/det(mat);
      temp_c=-det(mat3)/det(mat);
      temp_r=det(mat4)/det(mat);

      a=temp_a/-2;
      b=temp_b/-2;
      c=temp_c/-2;
      r=sqrt(a^2+b^2+c^2-temp_r);
      
      xx=[]; 
      nCurInlyerCount=0;        %初始化点数为0个
     %%  是否符合模型？ 
     for k=1:nDataLen
        CurModel=[a b c r];
%         pdist=sqrt((Data(1,k)-centerx).^2+(Data(2,k)-centery).^2);
        pdist=sqrt((Data(1,k)-a).^2+(Data(2,k)-b).^2+(Data(3,k)-c).^2);
        CurMask =(abs(r-pdist)< dThreshold);     %到直线距离小于阈值的点符合模型,标记为1
        nCurInlyerCount =nCurInlyerCount+CurMask;             %计算符合椭圆模型的点的个数
        if(CurMask==1)
             xx =[xx,Data(:,k)];
        end
     end
       %% 选取最佳模型
        if nCurInlyerCount > nMaxInlyerCount   %符合模型的点数最多的模型即为最佳模型
            nMaxInlyerCount = nCurInlyerCount;
            Ellipse_mask = CurMask;
             Ellipse_model = CurModel;
             Ellipse_points = [A B C D];
             Ellipse_x =xx;
        end
     
end
p=Ellipse_model;
% %% 由符合点拟合圆
% %椭圆一般方程：Ax2+Bxy+Cy2+Dx+Ey+F=0
% x=Ellipse_x';
% F=@(p)(x(:,1)-p(1)).^2 + (x(:,2)-p(2)).^2 + (x(:,3)-p(3)).^2- ball_radius^2;
% % F=@(p)(x(:,1)-p(1)).^2 + (x(:,2)-p(2)).^2 - p(3)^2;
% % p0=[1,1,1];
% p=lsqnonlin(F,[ball(1),ball(2),ball(3)]');%拟合的参数