close all;clear;%clc;
%% ��ȡ����
algorithem1 ="SFKNN��outlier�ķ���";

dataname="jain.txt";KK=2;

path=fullfile('D:\���ݼ�\',dataname);
AA = load(path);
%AA=load('D:\���ݼ�\Wine.txt');KK=3;

%% ��ѽ����¼
%��ѽ����¼��markdown�ļ���"D:\Markdown file\������ѽ����¼.md"

%% �����й�һ�� 
a=find(AA(:,end)==0);
AA(a,end)=3;
data=AA(:,1:end-1);
data= libsvmscale(data,0,1);  % ���ݹ�һ��
[rows,dim]=size(AA);
A=[data,AA(:,end)];
%% ����ʵͼ��
ShowClusterA(A,'origin graphic')
answer=AA(:,end);  %��ʵ��ǩ
%�������û�ã���Ϊ���㷨�������ñ�׼���Ȩ�������ˡ�
dist1=pdist2(data,data,'minkowski',2); %�ɿɷ�˹�������������ŷ����þ���
N=size(data,1);
%dist1=pdist(data);

K=2:50;
%K=15;
for p=1:length(K)  %KΪ���ڸ���
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

%% ��ͼ
% F=[B,cl'];
%if( dim==3)
ShowClusterA(F,'������');   
% ����������
scatter(data(bestresult.center,1),data(bestresult.center,2),100,'kh','MarkerFaceColor','w');%�ڱ�������


dashboard=table([result(:).ami]',[result(:).ari]',[result(:).fmi]',K','VariableNames',{'AMI','ARI','FMI','K'});
resultBestari = dashboard(dashboard.ARI==max([result(:).ari]),:)
