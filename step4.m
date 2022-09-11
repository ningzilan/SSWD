
clear;
clc;
tab=readtable('..\code\Biase_subcluster.csv');
CVI='Sil';
load("Biase.mat");
demension=1;
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
        else 
         tm=cumsum(latent)./sum(latent);
         ind=length(tm(tm<Thr));
        end
        temp=table(name,K_means_3(score(:,1:ind)));   
    end
      sub_cluster{i}=table2cell(temp); 
    
 end
 

 [nrow,ncol]=size(data);
 

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

Mat=zeros(nrow,nrow);
mm=length(sub_biny);
for i=1:mm
    if ~isempty(sub_biny{1,i})
    Mat=Mat+sub_biny{1,i};
    end
end

eva1=cell(10,1); OptionlK1=zeros(10,1);
  for j=1:10
    eva1{j,1}=evalclusters(Mat,'kmeans',CVI,'KList',1:10);
    OptionlK1(j,1)=eva1{j, 1}.OptimalK;
   end
   indx=tabulate(OptionlK1);
   [~,K_optimal]=max(indx(:,2))
   Y_optimal=kmedoids(Mat,K_optimal);
   AR = RandIndex(y,Y_optimal)
   NMI = NMI_AMI(y', Y_optimal')


    
       
 