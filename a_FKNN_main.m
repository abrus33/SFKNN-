close all;clear;%clc;
%% 读取数据
algorithem1 ="SFKNN改outlier改分配";

dataname="jain.txt";KK=2;

path=fullfile('D:\数据集\',dataname);
AA = load(path);
%AA=load('D:\数据集\Wine.txt');KK=3;

%% 最佳结果记录
%最佳结果记录在markdown文件中"D:\Markdown file\聚类最佳结果记录.md"

%% 数据列归一化 
a=find(AA(:,end)==0);
AA(a,end)=3;
data=AA(:,1:end-1);
data= libsvmscale(data,0,1);  % 数据归一化
[rows,dim]=size(AA);
A=[data,AA(:,end)];
%% 画真实图像
ShowClusterA(A,'origin graphic')
answer=AA(:,end);  %真实标签
%这里距离没用，因为在算法中重新用标准差加权距离算了。
dist1=pdist2(data,data,'minkowski',2); %闵可夫斯基的特殊情况：欧几里得距离
N=size(data,1);
%dist1=pdist(data);

K=2:50;
%K=15;
for p=1:length(K)  %K为近邻个数
  %If you want to choose centers manually, set AutoPick to 0, otherwise, number of centers
  for i=1
  tic
  result(p)=KnnDpc(data,answer,K(p),'AutoPick',KK,'Distance',dist1,'Ui',true);
  toc
  t(i)=toc;
  end

end
all_ari=[result.ari];
max_ = (all_ari==max(all_ari));
bestresult=result(max_);

all_K=[bestresult.K];
max_K=(all_K==max(all_K));
bestresult=bestresult(max_K);
F=[data,bestresult.cluster];
% [RI,ARI,NMI]=evolution(A,F,KK);
% acc=ClusteringAccuracy(result.cluster',AA(:,end));
% acc
% RI
% B=AA(:,1:2);
% F=[data,result.cluster];
% if( dim==3)
% ShowClusterA(F,'DPC-CNN-MCM on Pathbased');
% hold on
% scatter(data(result.center,1),data(result.center,2),'h','k');
% end

%% 画图
% F=[B,cl'];
%if( dim==3)
ShowClusterA(F,'聚类结果');   
% 画聚类中心
scatter(data(bestresult.center,1),data(bestresult.center,2),100,'kh','MarkerFaceColor','w');%黑边六角星


dashboard=table([result(:).ami]',[result(:).ari]',[result(:).fmi]',K','VariableNames',{'AMI','ARI','FMI','K'});
resultBestari = dashboard(dashboard.ARI==max([result(:).ari]),:)
