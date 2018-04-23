%% before running this: load desired posterior variables


%%
% JND differences between group and subject models FOR SUBJECT 1
tmp_FP1 = JND_group(:,1,1) - JND_subj(:,1);
% hist(tmp_FP1);
tmp_FP2 = JND_group(:,1,2) - JND_subj(:,1);
% hist(tmp_FP2);
tmp_FP3 = JND_group(:,1,3) - JND_subj(:,1);
% hist(tmp_FP3);
tmp_FP4 = JND_group(:,1,4) - JND_subj(:,1);
% hist(tmp_FP4);

TMP =zeros(20000,4);
TMP(:,1) = tmp_FP1;
TMP(:,2) = tmp_FP2;
TMP(:,3) = tmp_FP3;
TMP(:,4) = tmp_FP4;

for i=1:4
plotPosteriorDensity(TMP(:,i),'JND Compare Models', [2,2,i], 1, []);
end;
JNDpost_diff = TMP;

%%
% PSE differences between group and subject models FOR SUBJECT 1
tmp_FP1 = PSE_group(:,1,1) - PSE_subj(:,1);
% hist(tmp_FP1);
tmp_FP2 = PSE_group(:,1,2) - PSE_subj(:,1);
% hist(tmp_FP2);
tmp_FP3 = PSE_group(:,1,3) - PSE_subj(:,1);
% hist(tmp_FP3);
tmp_FP4 = PSE_group(:,1,4) - PSE_subj(:,1);
% hist(tmp_FP4);

TMP =zeros(20000,4);
TMP(:,1) = tmp_FP1;
TMP(:,2) = tmp_FP2;
TMP(:,3) = tmp_FP3;
TMP(:,4) = tmp_FP4;

for i=1:4
plotPosteriorDensity(TMP(:,i),'PSE Compare Models', [2,2,i], 1, []);
end;
PSEpost_diff = TMP;

save('PSEpost_diff.mat','PSEpost_diff');

%% Compare conditions FOR SUBJECT 1 (for now)

% choose one of these...
variablename = JND_group;
% variablename = JND_subj;
% variablename = PSE_group;
% variablename = PSE_subj;

% JND differences between groups 
tmp_FP12 = JND_group(:,1,1) - JND_group(:,1,2);
tmp_FP13 = JND_group(:,1,1) - JND_group(:,1,3);
tmp_FP14 = JND_group(:,1,1) - JND_group(:,1,4);
tmp_FP23 = JND_group(:,1,2) - JND_group(:,1,3);
tmp_FP24 = JND_group(:,1,2) - JND_group(:,1,4);
tmp_FP34 = JND_group(:,1,3) - JND_group(:,1,4);

TMP =zeros(20000,6);
TMP(:,1) = tmp_FP12;
TMP(:,2) = tmp_FP13;
TMP(:,3) = tmp_FP14;
TMP(:,4) = tmp_FP23;
TMP(:,5) = tmp_FP24;
TMP(:,6) = tmp_FP34;

for i=1:6
plotPosteriorDensity(TMP(:,i),'Compare Conditions', [2,3,i], 1, []);
end;
JNDpost_GRPdiff = TMP;
save('JNDpost_GRPdiff_S1.mat','JNDpost_GRPdiff');


%%
% choose one of these...
% variablename = JND_group;
% variablename = JND_subj;
variablename = PSE_group;
% variablename = PSE_subj;

% PSE differences between groups 
tmp_FP12 = PSE_group(:,1,1) - PSE_group(:,1,2);
tmp_FP13 = PSE_group(:,1,1) - PSE_group(:,1,3);
tmp_FP14 = PSE_group(:,1,1) - PSE_group(:,1,4);
tmp_FP23 = PSE_group(:,1,2) - PSE_group(:,1,3);
tmp_FP24 = PSE_group(:,1,2) - PSE_group(:,1,4);
tmp_FP34 = PSE_group(:,1,3) - PSE_group(:,1,4);

TMP =zeros(20000,6);
TMP(:,1) = tmp_FP12;
TMP(:,2) = tmp_FP13;
TMP(:,3) = tmp_FP14;
TMP(:,4) = tmp_FP23;
TMP(:,5) = tmp_FP24;
TMP(:,6) = tmp_FP34;

for i=1:6
plotPosteriorDensity(TMP(:,i),'PSE Compare Conditions', [2,3,i], 1, []);
end;
PSEpost_GRPdiff = TMP;
save('PSEpost_GRPdiff_S1.mat','PSEpost_GRPdiff');
