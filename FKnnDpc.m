function [result]=FKnnDpc(data,answer,K,varargin)
%
%% ������ʼ��
%Parse arguments 
parser=inputParser;                %����һ�� inputParser �������ڽ����������
parser.CaseSensitive=false;        %��Сд������
parser.KeepUnmatched=true;         % ����δƥ��Ĳ���
parser.StructExpand=true;
addRequired(parser,'data',@isreal);
addRequired(parser,'answer',@(answer)length(answer)==size(data,1));
addRequired(parser,'K',@(K)K>=0&&K<=length(answer));
addParameter(parser,'AutoPick',0,@(center)center>=0&&center<=length(answer));
addParameter(parser,'Distance',[],@(dist1)all(~diag(dist1))&&issymmetric(dist1));
addParameter(parser,'Ui',false);
parse(parser,data,answer,K,varargin{:}); %������ɺ󣬲�����ֵ������result�ṹ����

import java.util.LinkedList
import Library.*

N=size(data,1);
ns=N;
%% ���ݹ�һ��
%Normalization
% data=(data-min(data))./(max(data)-min(data));
% data(isnan(data))=0;
data= libsvmscale(data,0,1);  % ���ݹ�һ��
%Calculate dist1
 if isempty(parser.Results.Distance)%�ж��Ƿ�Ϊ��
  dist1=squareform(pdist(data)); %pdist��ÿ�����ݵ�֮��ľ������Ϊһ�У�squareform���Խ�pdist��������������ɾ��󣬵ڼ��еڼ������ö�Ӧ���������ݵ�֮��ľ���
  %%�ñ�׼����빫ʽ
%   s = std(data,1,1);
%   w = s/sum(s);
%   weighted_data = data.*w;
%   dist1 = squareform(pdist(weighted_data));
 else
   dist1=parser.Results.Distance;
 end

%% ����MNN
[knn_dist,knn_x]=sort(dist1,2);
knn_x=knn_x(:,2:K+1);
knn_dist=knn_dist(:,2:K+1);
adj = zeros(N,N);
for i=1:N
    adj(i,knn_x(i,:)) = 1; %���i��k����
end

mutual_adj = adj & adj' ; %���㻥���ھ���
%ת��Ϊ��Ԫ����洢���
%mutual_neighbors = arrayfun(@(x) find(mutual_adj(x,:)), 1:ns, 'UniformOutput', false);
%%ת��Ϊ�д洢���� =====
% Ԥ�����ڴ棨zero��䣩
mutual_counts = sum(mutual_adj, 2);       % ÿ����Ļ�������
max_mutual = max(mutual_counts);          % ��󻥽�����
mnn_matrix = NaN(N, max_mutual);       % ��ʼ���������
mnn_dist=NaN(N,max_mutual); %��ʼ��Mnn�������
% �������
for i = 1:ns
    neighbors = find(mutual_adj(i,:));    % �ҳ���ǰ��Ļ�����
    mnn_matrix(i, 1:length(neighbors)) = neighbors;
    mnn_dist(i,1:length(neighbors)) = dist(i,neighbors);
end

%% ��rho
%Calculate \rho

rho=sum(exp(-knn_dist),2);    

%% ��delta
%Calculate \delta
delta=inf(1,N);
deltaSelect=zeros(1,N);
[~,rhoOrder]=sort(rho,'descend');    %�ܶȽ���
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

%% ��gamma
%Calculate \gamma
gamma=rho.*delta';

%% ѡ���ܶ����ĵ㣬�ó���Ϊ�мල�����ĵ�
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

%% ����Ⱥ��
%Find outlier
theta=mean(knn_dist(:,K)); %ÿ���㵽k���ڵ��������ƽ��ֵ
%0 means outlier%���Ϊ0�ĵ�Ϊ��Ⱥ��
cluster(knn_dist(:,K)>theta)=0;
cluster(center)=1:NC;%�ܶȷ�ֵ����Ϊ1��NC������Ⱥ����λ-1����Ⱥ����Ϊ0

% [~,N] = size(find(cluster==0)) ;
% F=[data,cluster'];
% ShowClusterA(F,'��Ⱥ��');

%% ���ĵ����
%Assign strategy 1
 queue=LinkedList; 
 for p=center
   %һ�磬�ܲ���
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
% ShowClusterA(F,'һ�׶η����');


%% ������Ⱥ��
%% ������Ⱥ���δ����ĺ��ĵ�
%��ԭ�㷨�����𣺼���������ʱֻ�Ǽ�ͳ��k'�������ж�������ĳ���أ�
%����������ʱ��Ҳֻ������ͳ�ƣ�ûʲô�ر�Ĺ�ʽ
unas=find(cluster<=0); %δ�����
unasCount=length(unas);
while true
  recog=zeros(unasCount,NC);
  for p=1:unasCount  %��������δ����ĵ�p
    for o=1:K    %��������δ������k����
      u=cluster(knn_x(unas(p),o)); %��k�������Ƿ��б������
      if u>0
        recog(p,u)=recog(p,u)+1; %����У�p������ص������Ⱦͼ�1
      end
    end
  end
  whichK=max(recog(:));
  if whichK==K %���ĳ���������k��������ͬһ����
    %[whichPoint,whichCluster]=ind2sub(size(recog),find(recog==whichK));
    [whichPoint,whichCluster]=find(recog==whichK);
     cluster(unas(whichPoint))=whichCluster;
    unas=find(cluster<=0);
    unasCount=length(unas);
  elseif whichK>0 %���ĳ������k���ڱ�����
    %[whichPoint,whichCluster]=ind2sub(size(recog),randsample(find(recog==whichK),1));
    [whichPoint,whichCluster]=find(recog==whichK);
    cluster(unas(whichPoint(1,1)))=whichCluster(1,1);
    unas=find(cluster<=0);
    unasCount=length(unas);
  else
    break;
  end
end
%% ���������
%Assign outlier
%����������䵽���������k���ڵĴ���
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
%�������δ����ĵ㣬��δ����ĵ���Ϊһ���´�
if min(cluster)<=0
    cluster(cluster <=0) = max(cluster)+1;
end

%% ����
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