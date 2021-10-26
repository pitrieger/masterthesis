function [ B, F, loss ] = clusterwiseSCAPwithpartition( Zsup, I, K, Q, C,
partvec )
% Calculate the SCA-P solutions for each cluster with fixed partition
% I - number of subjects (Blocks)
% K - vector with number of observations for subjects in it
% Q - number of components
% C - number of clusters
% partvec - partition vector
kcum=zeros(I,2);
for n=1:I,
  kcum(n,1)=(sum(K(1:n-1))+1);
  kcum(n,2)=(sum(K(1:n)));
end;

B=[ ];
F=zeros(sum(K),C*Q); % Initialization of Fsup
for c=1:C % Perform SCA-P per cluster
  subj_c=find(partvec==c)';
  K_c=K(subj_c); % A vector with the K_n's for the subjects assigned to cluster cl
  zsup_c=[ ];
  kcum_c=zeros(length (subj_c),2);
  for n_c=1:length(subj_c)
    Z=Zsup(kcum(subj_c(n_c),1):kcum(subj_c(n_c),2),:);
    Zsup_c=[Zsup_c;Z];
    kcum_c(n_c,1)=(sum(K_c(1:n_c-1))+1);
    kcum_c(n_c,2)=(sum(K_c(1:n_c)));
  end;

  [U,S,V]=svd(Zsup_c,0);
  cB=sqrt(1/sum(K_c))*V(:,1:Q)*S(1:Q,1:Q);
  B=[B,cB];
  for n_c=1:length(subj_c)
    F(kcum(subj_c(n_c),1):kcum(subj_c(n_c),2),(c-1)*Q+1:c*Q)=sqrt(sum(K_c))*U(kcum_c(n_c,1):kcum_c(n_c,2),1:Q);
  end
end
% Find the loss of a model
FB = F*B';
loss = ssq(Zsup - FB);
end