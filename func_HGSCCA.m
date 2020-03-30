function  [accuracy,spec,sens,PPV,NPV,AUC]=func_HGSCCA(train_epi,test_epi,train_miRNA,test_miRNA,train_RNA,test_RNA,train_label,test_label,opts,k)
alpha=opts.alpha;
r_e=opts.r_e;
r_mi=opts.r_mi;
r_R=opts.r_R;
beta_e=opts.beta_e;
beta_mi=opts.beta_mi;
beta_R=opts.beta_R;
currentFolder = pwd;
addpath(genpath(currentFolder))
eps=1e-10;
[H_epi,Lh_epi] = cons_hypergraph(train_epi, k); 
[H_miRNA,Lh_miRNA] = cons_hypergraph(train_miRNA, k); 
[H_RNA,Lh_RNA] = cons_hypergraph(train_RNA, k); 

train_data=[train_epi,train_miRNA,train_RNA];
test_data=[test_epi,test_miRNA,test_RNA];

nTrain=size(train_data,1);
nTest=size(test_data,1);

w_e=ones(size(train_epi,2),1);
w_mi=ones(size(train_miRNA,2),1);
w_R=ones(size(train_RNA,2),1);
G=zeros(size(train_label,1),1);



%%%%%%%%%% update G %%%%%%%%%
for i=1:10
  G=alpha*1/3*(train_epi*w_e+train_miRNA*w_mi+train_RNA*w_R);

%%%%%%%%%%%update w_e, w_mi,w_RNA%%%%%%%%%%%%%%%%%%
  D_e=diag(1./sqrt(w_e.*w_e+1e-4));
  D_mi=diag(1./sqrt(w_mi.*w_mi+eps));
  D_R=diag(1./sqrt(w_R.*w_R+1e-4));

  w_e=((alpha+1/nTrain)*train_epi'*train_epi+beta_e*train_epi'*Lh_epi*train_epi+r_e*D_e)\(alpha*train_epi'*G+(1/nTrain)*train_epi'*train_label);
  w_mi=inv((alpha+1/nTrain)*train_miRNA'*train_miRNA+beta_mi*train_miRNA'*Lh_miRNA*train_miRNA+r_mi*D_mi)*(alpha*train_miRNA'*G+(1/nTrain)*train_miRNA'*train_label);
  w_R=inv((alpha+1/nTrain)*train_RNA'*train_RNA+beta_R*train_RNA'*Lh_RNA*train_RNA+r_R*D_R)*(alpha*train_RNA'*G+(1/nTrain)*train_RNA'*train_label);
  ZR(i,:)=w_R';
  ZE(i,:)=w_e';
  ZMi(i,:)=w_mi';
  f_val(i)=alpha*(norm(G-train_epi*w_e,2)+norm(G-train_miRNA*w_mi,2)+norm(G-train_RNA*w_R,2))+1/nTrain*(norm(train_label-train_epi*w_e,2)+norm(train_label-train_miRNA*w_mi,2)+norm(train_label-train_RNA*w_R,2))+r_e*norm(w_e,1)+r_mi*norm(w_mi,1)+r_R*norm(w_R,1)+...
      beta_e*w_e'*train_epi'*Lh_epi*train_epi*w_e+ beta_mi*w_mi'*train_miRNA'*Lh_miRNA*train_miRNA*w_mi+beta_R*w_R'*train_RNA'*Lh_RNA*train_RNA*w_R;
  
end
ind_Epi=find(abs(w_e)>=1e-3);
ind_RNA=find(abs(w_R)>=1e-3);
ind_miRNA=find(abs(w_mi)>=1e-3);

len_E=length(ind_Epi);
len_R=length(ind_RNA);
len_M=length(ind_miRNA);

s_train_Epi=train_epi(:,ind_Epi);
s_train_RNA=train_RNA(:,ind_RNA);
s_train_miRNA=train_miRNA(:,ind_miRNA);

s_test_Epi=test_epi(:,ind_Epi);
s_test_RNA=test_RNA(:,ind_RNA);
s_test_miRNA=test_miRNA(:,ind_miRNA);

s_train=[s_train_Epi,s_train_RNA,s_train_miRNA];
s_test=[s_test_Epi,s_test_RNA,s_test_miRNA];

parameter=['-t 2',' -c ',num2str(4), ' -g ',num2str(2^-10)];
model=svmtrain(double(train_label),double(s_train),parameter);
[predict_b, tempAccuracy, dec]=svmpredict(double(test_label), double(s_test), model);
[X1,Y1,T,AUC] = perfcurve(test_label,dec,1);
[accuracy,spec,sens,PPV,NPV]=getClassificationResult(predict_b,test_label);

  
end 











