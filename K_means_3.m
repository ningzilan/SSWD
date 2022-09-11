function [Y_optimal,K_optimal,W_optimal,RCH_max,RCH_K] = K_means_3(data)
[nrow,mcol]=size(data);      
Sdata=standard_z(data);   
[E0,~]=similarity_euclid(Sdata);
E1=nornaliz(E0);
R0=corrcoef(Sdata');
R1=nornaliz(R0);
 for K=2:10
	[RCH_K_Y(:,K),RCH_K(1,K),w_optimal(1,K)]=specify_K_means3(Sdata,E1,R1,K);
  end
[RCH_max,lab]=max(RCH_K);
K_optimal=lab;
Y_optimal=RCH_K_Y(:,lab);
W_optimal=w_optimal(1,lab);
end
	
function [sdata]=standard_z(data)
	nrow=size(data,1);
	colmean=mean(data);
	data=data-repmat(colmean,nrow,1);
	sdata=data./std(data);
end

function [Edist,dmax,dmin]=similarity_euclid(data)
	nrow=size(data,1);
	Edist=zeros(nrow,nrow);
	dmax=0;
	dmin=0;
	for i=1:nrow
		for j=1:nrow
			Edist(i,j)=norm(data(i,:)-data(j,:));
		end
	end
end

function [M]=nornaliz(data)
[nrow,mcol]=size(data);
M=data;
M(logical(eye(size(M))))=0;
dmax=max(max(M));
M(logical(eye(size(M))))=5;
dmin=min(min(M));
cmax=dmax-dmin;
M = M - repmat(dmin,nrow,1);
M = M./repmat(cmax,nrow,1);
M(logical(eye(size(M))))=0;
end

function [RCH_K_Y,RCH_K,w_optimal]=specify_K_means3(Sdata,E,R,K)
[nrow,mcol]=size(Sdata);
C=zeros(K,mcol);
CH_W_Y=zeros(nrow,11);
CH=zeros(1,11);
n_w=1;
for w=0:0.1:1
	ED_Pearson=w*E+(1-w)*(1-R);
	ED_Pearson(logical(eye(size(ED_Pearson))))=0;
	[c1,c2]=find(ED_Pearson==max(max(ED_Pearson)));
    center=[c1(1,1),c2(1,1)];
    for v=3:K
		dis0=ED_Pearson(center,:);
		[d3,~]=min(dis0);
		[~,d3]=max(d3);
		center(v)=d3;
    end
	dis=ED_Pearson(center,:);
	[~,fy]=min(dis);
	count=0;
	C=Sdata(center,:);
	C_adjust=zeros(K,mcol);
    center_adjust=zeros(1,K);
	while(count<50 && ~isequal(C,C_adjust))
		C_adjust=C;
        center_adjust=center;
		for i=1:K
			Cindex=find(fy==i);
			n_i=length(Cindex);
			if n_i>1        
				data_i=Sdata(Cindex,:);
				ED_i=ED_Pearson(Cindex,:);
				ED_i=ED_i(:,Cindex);
				[~,lab]=min(sum(ED_i,2));
				C(i,:)=data_i(lab,:);
                center(1,i)=lab;
			end
        end
        fy=clustering(Sdata,C,w);  
        count=count+1;
    end
    [CH_W(n_w),CH_W_0(n_w)]=get_CH(Sdata,fy);
    CH_W_Y(:,n_w)=fy;
    n_w=n_w+1;
end

[CH_Wmax,lab]=max(CH_W);
RCH_K=CH_Wmax/finv(0.95,K-1,nrow-K);
RCH_K_Y=CH_W_Y(:,lab);
w_optimal=(lab-1)*0.1;
end

function [fy]=clustering(Sdata,C,w)
[nrow,mcol]=size(Sdata);
[K,~]=size(C);
Sdata=Sdata';
C=C';
	Edist=size(nrow,K);
	for j=1:K
		for i=1:nrow
			Edist_c(i,j)=norm(Sdata(:,i)-C(:,j));
			r=corrcoef(Sdata(:,i),C(:,j));
            R_c(i,j)=r(1,2);
		end
	end
	[Edist_c0]=nornaliz_1(Edist_c);
	[R_c0]=nornaliz_1(R_c);
	ED_Pearson_1=w*Edist_c0+(1-w)*(1-R_c0);
	for i=1:nrow
	[~,lab]=min(ED_Pearson_1(i,:));
	fy(i,1)=lab;
    end
end


function  [M]=nornaliz_1(D)
[nrow,mcol]=size(D);
dmax=max(max(D));
dmin=min(min(D));
cmax=dmax-dmin;
M=D-repmat(dmin,nrow,1);
M=M./repmat(cmax,nrow,1);
end

function [CH,CH_0]=get_CH(Sdata,fy)
[nrow,mcol]=size(Sdata);
for i=1:mcol
	[~,table,~]=anova1(Sdata(:,i),fy,'off');
    F(i)=table{2,5};
	Mst(i)=table{2,4};
	Mse(i)=table{3,4};
	dft=table{2,3};
	dfe=table{3,3};
end
CH=sum(Mst)/sum(Mse);
CH_0=sum(F);
end
	