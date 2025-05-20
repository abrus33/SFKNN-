function [result]=FKnnDpc(data,answer,K,varargin)
%
%% 参数初始化
%Parse arguments 
parser=inputParser;                %创建一个 inputParser 对象，用于解析输入参数
parser.CaseSensitive=false;        %大小写不敏感
parser.KeepUnmatched=true;         % 保留未匹配的参数
parser.StructExpand=true;
addRequired(parser,'data',@isreal);
addRequired(parser,'answer',@(answer)length(answer)==size(data,1));
addRequired(parser,'K',@(K)K>=0&&K<=length(answer));
addParameter(parser,'AutoPick',0,@(center)center>=0&&center<=length(answer));
addParameter(parser,'Distance',[],@(dist1)all(~diag(dist1))&&issymmetric(dist1));
addParameter(parser,'Ui',false);
parse(parser,data,answer,K,varargin{:}); %解析完成后，参数的值保存在result结构体中

import java.util.LinkedList
import Library.*

N=size(data,1);
ns=N;
%% 数据归一化
%Normalization
% data=(data-min(data))./(max(data)-min(data));
% data(isnan(data))=0;
data= libsvmscale(data,0,1);  % 数据归一化
%Calculate dist1
 if isempty(parser.Results.Distance)%判断是否为空
  dist1=squareform(pdist(data)); %pdist将每个数据点之间的距离计算为一行，squareform可以将pdist出来的行向量变成矩阵，第几行第几列正好对应哪两个数据点之间的距离
  %%用标准差距离公式
%   s = std(data,1,1);
%   w = s/sum(s);
%   weighted_data = data.*w;
%   dist1 = squareform(pdist(weighted_data));
 else
   dist1=parser.Results.Distance;
 end

%% 改用MNN
[knn_dist,knn_x]=sort(dist1,2);
knn_x=knn_x(:,2:K+1);
knn_dist=knn_dist(:,2:K+1);
adj = zeros(N,N);
for i=1:N
    adj(i,knn_x(i,:)) = 1; %标记i的k近邻
end

mutual_adj = adj & adj' ; %计算互近邻矩阵
%转换为单元数组存储结果
%mutual_neighbors = arrayfun(@(x) find(mutual_adj(x,:)), 1:ns, 'UniformOutput', false);
%%转换为行存储矩阵 =====
% 预分配内存（zero填充）
mutual_counts = sum(mutual_adj, 2);       % 每个点的互近邻数
max_mutual = max(mutual_counts);          % 最大互近邻数
mnn_matrix = NaN(N, max_mutual);       % 初始化结果矩阵
mnn_dist=NaN(N,max_mutual); %初始化Mnn距离矩阵
% 填充数据
for i = 1:ns
    neighbors = find(mutual_adj(i,:));    % 找出当前点的互近邻
    mnn_matrix(i, 1:length(neighbors)) = neighbors;
    mnn_dist(i,1:length(neighbors)) = dist(i,neighbors);
end

%% 求rho
%Calculate \rho

rho=sum(exp(-knn_dist),2);    

%% 求delta
%Calculate \delta
delta=inf(1,N);
deltaSelect=zeros(1,N);
[~,rhoOrder]=sort(rho,'descend');    %密度降序
for p=2:N
  for o=1:p-1
    if dist1(rhoOrder(p),rhoOrder(o))<delta(rhoOrder(p))
      delta(rhoOrder(p))=dist1(rhoOrder(p),rhoOrder(o));
      deltaSelect(rhoOrder(p))=rhoOrder(o);
    end
  end
end
delta(rhoOrder(1))=-1;
delta(rhoOrder(1))=max(delta);

%% 求gamma
%Calculate \gamma
gamma=rho.*delta';

%% 选择密度中心点，该程序为有监督求中心点
%Ui
%If AutoPick is set, I assume you want to make a benchmark and will decide centers automatically
if parser.Results.AutoPick
  gammaSort=sort(gamma,'descend');
  x=nan;
  y=mean(gammaSort(parser.Results.AutoPick:parser.Results.AutoPick+1));
  click=2;
else
  [x,y,click,fig,ax1,ax2]=UiMakeDecision(rho,delta);
end

%Assign center
%Cluster to which it belongs, notate unassigned in -1
cluster=-ones(1,N)';
%Id of centers
center=[];
if click==1
  %x->rhoLeast, y->deltaLeast
  center=intersect(find(rho>x),find(delta>y));
elseif click==2
  %y->gammaLeast
  center=find(gamma>y);
end
%Number of clusters
NC=length(center);

%% 找离群点
%Find outlier
theta=mean(knn_dist(:,K)); %每个点到k近邻的最大距离的平均值
%0 means outlier%标记为0的点为离群点
cluster(knn_dist(:,K)>theta)=0;
cluster(center)=1:NC;%密度峰值点标记为1：NC，非离群点标记位-1，离群点标记为0

% [~,N] = size(find(cluster==0)) ;
% F=[data,cluster'];
% ShowClusterA(F,'离群点');

%% 核心点分配
%Assign strategy 1
 queue=LinkedList; 
 for p=center
   %一坨，跑不了
   for i=1:K
   cluster(knn_x(p,i))=cluster(p);
   end
   for o=knn_x(p,1:K)
      queue.offerFirst(o);
   end
   while ~queue.isEmpty()
     this=queue.pollLast();
     for i=1:length(this)
      for j=1:K
     for next=knn_x(this(i),j)
         
       if cluster(next)<0 && dist1(this(i),next)<=mean(knn_dist(next,:))
         cluster(next)=cluster(this(i));
         queue.offerFirst(next);
       end
     end
     end
     end
   end
  end

%% Assign strategy 1
%  queue=LinkedList;
%  for p=1:length(center)   
%       queue.offerFirst(center(p)); 
%  end
% 
%    while ~queue.isEmpty()
%      this=queue.pollLast();
%      for next=knn_x(this,1:K)
%        if cluster(next)<0&&dist1(this,next)<=mean(knn_dist(this,1:K))
%          cluster(next)=cluster(this);
%          queue.offerFirst(next);
%        end
%      end
%    end
 
% [~,NN] = size(find(cluster==0)) ;
% F=[data,cluster'];
% ShowClusterA(F,'一阶段分配后');


%% 分配离群点
%% 分配离群点和未分配的核心点
%与原算法的区别：计算隶属度时只是简单统计k'近邻中有多少属于某个簇，
%更新隶属度时，也只是重新统计，没什么特别的公式
unas=find(cluster<=0); %未分配的
unasCount=length(unas);
while true
  recog=zeros(unasCount,NC);
  for p=1:unasCount  %遍历所有未分配的点p
    for o=1:K    %遍历所有未分配点的k近邻
      u=cluster(knn_x(unas(p),o)); %找k近邻中是否有被分配的
      if u>0
        recog(p,u)=recog(p,u)+1; %如果有，p对这个簇的隶属度就加1
      end
    end
  end
  whichK=max(recog(:));
  if whichK==K %如果某个点的所有k近邻属于同一个簇
    %[whichPoint,whichCluster]=ind2sub(size(recog),find(recog==whichK));
    [whichPoint,whichCluster]=find(recog==whichK);
     cluster(unas(whichPoint))=whichCluster;
    unas=find(cluster<=0);
    unasCount=length(unas);
  elseif whichK>0 %如果某个点有k近邻被分配
    %[whichPoint,whichCluster]=ind2sub(size(recog),randsample(find(recog==whichK),1));
    [whichPoint,whichCluster]=find(recog==whichK);
    cluster(unas(whichPoint(1,1)))=whichCluster(1,1);
    unas=find(cluster<=0);
    unasCount=length(unas);
  else
    break;
  end
end
%% 噪声点分配
%Assign outlier
%将噪声点分配到离他最近的k近邻的簇中
unas=find(cluster<=0);
for p=1:length(unas)
  %1 is itself
  for o=knn_x(unas(p),:)
    if cluster(o)>0
      cluster(p)=cluster(o);
      break;
    end
  end
end

%Remain same shape with answer
center=center';
%如果存在未分配的点，将未分配的点视为一个新簇
if min(cluster)<=0
    cluster(cluster <=0) = max(cluster)+1;
end

%% 评价
%cluster(center)=1:NC;
% Evaluation
if answer
  ami=Ami(answer,cluster);
  ari=Library.GetAri(answer,cluster);
  fmi=Library.GetFmi(answer,cluster);
%   acc=ClusteringAccuracy(answer,cluster);
else
  ami=nan;
  ari=nan;
  fmi=nan;
%   acc=nan;
end

%Ui
% if parser.Results.Ui
%   if parser.Results.AutoPick
%     [fig,ax1,ax2]=UiShowDecision(rho,delta,x,y,click,K);
%   end
%   data=UiReduceDimension(data);
%   cmap=UiGetColormap(NC);
%   UiDrawMarker(NC,x,y,click,K,fig,ax1,ax2);
%   UiDrawCenter(center,rho,delta,cmap,fig,ax1,ax2);
%   UiPlotResultAndRho(cluster,center,rho,data,cmap,fig,subplot(2,2,3));
%   UiPlotResultAndDelta(cluster,center,delta,data,cmap,fig,subplot(2,2,4));
% end

%Wrap all data into a structure
result=struct;
result.answer=answer;
result.NC=NC;
result.K=K;
result.dist2=nan;
result.rho=rho;
result.delta=delta;
result.deltaSelect=deltaSelect;
result.gamma=gamma;
result.cluster=cluster;
result.center=center;
result.ami=ami;
result.ari=ari;
result.fmi=fmi;
% result.acc=acc;
result.x=x;
result.y=y;
result.click=click;

end