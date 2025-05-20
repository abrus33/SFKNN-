clear
data_dir1 = 'D:\数据集\UCI数据集\';
%function batch_clustering_test(data_dir, result_path, K_values)
    % 输入：
    %   data_dir    : 数据集目录（字符串，如 'datasets/'）
    %   result_path : 结果保存路径（字符串，如 'results/summary.xlsx'）
    %   K_values    : 近邻数K的候选值（数组，如 [5, 10, 15]）

    disp("SFKNN")
    disp("UCI数据集")

% 获取所有数据集文件
file_list = dir(fullfile(data_dir1,'*.txt'));

%遍历每个数据集
for fidx = 1:length(file_list)
    file_name = file_list(fidx).name;
    file_name_string=string(file_name);
    file_name_string=string(file_name);
    file_name_string=repmat(file_name_string,1,49);
    dataset_name = strrep(file_name,'.txt','');
    data_path = fullfile(data_dir1,file_name);
    AA = load(data_path);
    %输入真实簇类数，可以写一个数组记录每个数据集的簇类数，一起循环
   
    %% 数据归一化
    a=find(AA(:,end)==0);
    AA(a,end)=3;
    data = AA(:,1:end-1);
    data = libsvmscale(data,0,1);
    [rows,dim] = size(AA);
    A = [data,AA(:,end)];
    answer = AA(:,end);

    KK=max(answer);
    K=2:50;    
    for p = 1:length(K)
        %for i = 1
            %tic
            result(p)=SFKnnDpc(data,answer,K(p),'AutoPick',KK,'Ui',true);
            %toc
            %t(i) = toc;
        %end
  
    end
    all_ari = [result.ari];
    max_ = (all_ari==max(all_ari));
    bestresult = result(max_);

    all_K = [bestresult.K];
    max_K=(all_K==max(all_K));
    %bestresult=bestresult(K==max_K);
    dashboard=table(file_name_string',[result(:).ami]',[result(:).ari]',[result(:).fmi]',K','VariableNames',{'dataset','AMI','ARI','FMI','K'});
    resultBestari = dashboard(dashboard.ARI==max([result(:).ari]),:)
    
end

data_dir2 = 'D:\数据集\人工数据集\';
disp("人工数据集")
% 获取所有数据集文件
file_list = dir(fullfile(data_dir2,'*.txt'));

%遍历每个数据集
for fidx = 1:length(file_list)
    file_name = file_list(fidx).name;
    file_name_string=string(file_name);
    file_name_string=string(file_name);
    file_name_string=repmat(file_name_string,1,49);
    dataset_name = strrep(file_name,'.txt','');
    data_path = fullfile(data_dir2,file_name);
    AA = load(data_path);
    %输入真实簇类数，可以写一个数组记录每个数据集的簇类数，一起循环
   
    %% 数据归一化
    a=find(AA(:,end)==0);
    AA(a,end)=3;
    data = AA(:,1:end-1);
    data = libsvmscale(data,0,1);
    [rows,dim] = size(AA);
    A = [data,AA(:,end)];
    answer = AA(:,end);

    K=2:50;
    KK=max(answer);
    for p = 1:length(K)
        %for i = 1
            %tic
            result(p)=SFKnnDpc(data,answer,K(p),'AutoPick',KK,'Ui',true);
            %toc
            %t(i) = toc;
        %end
  
    end
    all_ari = [result.ari];
    max_ = (all_ari==max(all_ari));
    bestresult = result(max_);

    all_K = [bestresult.K];
    max_K=(all_K==max(all_K));
    %bestresult=bestresult(K==max_K);
    dashboard=table(file_name_string',[result(:).ami]',[result(:).ari]',[result(:).fmi]',K','VariableNames',{'dataset','AMI','ARI','FMI','K'});
    resultBestari = dashboard(dashboard.ARI==max([result(:).ari]),:)
    
end