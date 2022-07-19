%% SSWD  
%% Zilan Ning 
%% 2022/5/25
clear;
clc;
tab=readtable('C:\R\cellcluster\experiment_data_906\log2_10\k_2_10\sd1000\Li58_subcluster_sd.csv');
load("GSE_Li_58.mat");
CVI='Sil';
demension=1; 
k=length(unique(y));
cellname=tab{1:end-1,1};
data=tab(:,2:end);
data=table2array(data);
[nrow,~]=size(data);
cluster=max(data(nrow,:));
sub_cluster=cell(1,cluster);
 for i=1:cluster
    sub_data=data(1:end-1,data(end,:)==i);
    sub_all=table(cellname,sub_data);
    sub_all((find(sum(sub_data,2)==0)),:)=[];  
    name=sub_all{:,1}; 
    sub_data=table2array(sub_all(:,2:end));
    [subrow,subcol]=size(sub_data);
    if subrow<3 || subcol<3
          continue;
    else
        
        [coeff,score,latent,tsquare]=pca(zscore(sub_data));
        if demension==1 
        ind=Elbow(score);
        %[pc_data]=pc(sub_data);
        else 
         tm=cumsum(latent)./sum(latent);
         ind=length(tm(tm<Thr));
        end
        temp=table(name,K_means_3(score(:,1:ind)));   % 一致性聚类
        %temp=table(name,K_means_3(pc_data));
        %temp=score(:,1:ind);   %汇总聚类
    end
      sub_cluster{i}=table2cell(temp); %一致性聚类
     %sub_cluster{i}=temp;  %汇总聚类
 end
 
 %汇总聚类
%  Mat=[];
%  for i=1:cluster
%      Mat=[Mat,sub_cluster{i}];
%  end
%  %Mat=mapstd(Mat);
%  [Y_optimal,K_optimal] = K_means_3(Mat);
 
%一致性聚类
 [nrow,ncol]=size(data);
%全部子空间聚类
sub_biny=cell(1,cluster);
for i=1:cluster
    if isempty(sub_cluster{1,i})
        continue;
    else
    [subrow,subcol]=size(sub_cluster{1,i}(:,2));
     c1=cell2mat(sub_cluster{1,i}(:,2));
        for m=1:subrow
            for n=1:subrow
                if c1(m)==c1(n)
                    sub_biny{1,i}(n,m)=1;
                    sub_biny{1,i}(m,n)=1;
                end
            end
        end
    end
end

 %补齐缺失cell 
 [nrow,ncol]=size(data);
 nrow=nrow-1;    
 for i=1:cluster
    if isempty(sub_cluster{1,i})
         continue;
     else
       [subrow,subcol]=size(sub_cluster{1,i});
        if subrow<nrow
         sub_cellname=sub_cluster{1,i}(:,1);
         del=setdiff(cellname,sub_cellname);
         sub_cellname=[sub_cellname;del];
         [~,ind]=sort(sub_cellname);
         sdata=zeros(nrow,nrow);
         sdata(1:subrow,1:subrow)=sdata(1:subrow,1:subrow)+sub_biny{1,i};
         sdata=sdata(ind,:);
         sdata=sdata(:,ind);
         sub_biny{1,i}=sdata;
        end
     end
 end

%多个聚类结果根据CSPA规整为一个cell*cell矩阵
Mat=zeros(nrow,nrow);
mm=length(sub_biny);
for i=1:mm
    if ~isempty(sub_biny{1,i})
    Mat=Mat+sub_biny{1,i};
    end
end
% 
% %一致性聚类


  %PAM
      eva1=cell(10,1); OptionlK1=zeros(10,1);
      for j=1:10
         eva1{j,1}=evalclusters(Mat,'kmeans',CVI,'KList',1:10);
         OptionlK1(j,1)=eva1{j, 1}.OptimalK;
      end
       indx=tabulate(OptionlK1);
      [~,K_optimal]=max(indx(:,2))
       Y_optimal=kmedoids(Mat,K_optimal);
       AR = RandIndex(y,Y_optimal)
       [NMI,AMI] = NMI_AMI(y', Y_optimal')
       % save('Y_method3_Eblow_PAM_Biase.mat','Y_optimal');


    
       
 