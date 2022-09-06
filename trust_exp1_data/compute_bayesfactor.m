load 'modelfiiting_evidence.mat'

%% compute group bayes facor

hybrid_MF_GBF = exp(sum(T_Model_evidence_hybrid(:,1)) - sum(T_Model_evidence_MF_only(:,1)));
hybrid_MB_GBF = exp(sum(T_Model_evidence_hybrid(:,1)) - sum(T_Model_evidence_MB_only(:,1)));
MF_MB_GBF = exp(sum(T_Model_evidence_MF_only(:,1)) - sum(T_Model_evidence_MB_only(:,1)));

UN_hybrid_MF_GBF = exp(sum(Un_Model_evidence_hybrid(:,1)) - sum(Un_Model_evidence_MF_only(:,1)));
UN_hybrid_MB_GBF = exp(sum(Un_Model_evidence_hybrid(:,1)) - sum(Un_Model_evidence_MB_only(:,1)));
UN_MF_MB_GBF = exp(sum(Un_Model_evidence_MF_only(:,1)) - sum(Un_Model_evidence_MB_only(:,1)));



figure;
col = [0.7 0.7 0.7];
barh(T_Model_evidence_hybrid(:,1) - T_Model_evidence_MF_only(:,1), 'facecolor', col, 'edgecolor', 'none');
ylim([0 size(T_Model_evidence_hybrid,1)+1])
xlim([-5 10])
title('model evidence difference_honesty', 'FontSize', 16)
ylabel('subject', 'FontSize', 12)
xlabel('\Delta lme', 'FontSize', 12)
text(-5,15,'MF Model', 'FontSize', 10)
text(5,15,'Hybrid Model', 'FontSize', 10)
box off


figure;
col = [0.7 0.7 0.7];
barh(T_Model_evidence_hybrid(:,1) - T_Model_evidence_MB_only(:,1), 'facecolor', col, 'edgecolor', 'none');
ylim([0 size(T_Model_evidence_hybrid,1)+1])
xlim([-5 25])
title('model evidence difference_honesty', 'FontSize', 16)
ylabel('subject', 'FontSize', 12)
xlabel('\Delta lme', 'FontSize', 12)
text(-4,15,'MB Model', 'FontSize', 10)
text(15,15,'hybrid Model', 'FontSize', 10)
box off


% figure;
% col = [0.7 0.7 0.7];
% barh(T_Model_evidence_MF_only(:,1) - T_Model_evidence_MB_only(:,1), 'facecolor', col, 'edgecolor', 'none');
% ylim([0 size(T_Model_evidence_hybrid,1)+1])
% title('model evidence difference_honesty', 'FontSize', 16)
% ylabel('subject', 'FontSize', 12)
% xlabel('\Delta lme', 'FontSize', 12)
% text(-80,20,'MB', 'FontSize', 14)
% text(10,20,'MF', 'FontSize', 14)
% box off

figure;
col = [0.7 0.7 0.7];
barh(Un_Model_evidence_hybrid(:,1) - Un_Model_evidence_MF_only(:,1), 'facecolor', col, 'edgecolor', 'none');
ylim([0 size(T_Model_evidence_hybrid,1)+1])
xlim([-3 10])
title('model evidence difference_dishonesty', 'FontSize', 16)
ylabel('subject', 'FontSize', 12)
xlabel('\Delta lme', 'FontSize', 12)
text(-3,15,'MF Model', 'FontSize', 10)
text(7,15,'Hybrid Model', 'FontSize', 10)
box off


figure;
col = [0.7 0.7 0.7];
barh(Un_Model_evidence_hybrid(:,1) - Un_Model_evidence_MB_only(:,1), 'facecolor', col, 'edgecolor', 'none');
ylim([0 size(T_Model_evidence_hybrid,1)+1])
xlim([-3 15])
title('model evidence difference_dishonesty', 'FontSize', 16)
ylabel('subject', 'FontSize', 12)
xlabel('\Delta lme', 'FontSize', 12)
text(-3,15,'MB Model', 'FontSize', 10)
text(13,15,'hybrid Model', 'FontSize', 10)
box off

% 
% figure;
% col = [0.7 0.7 0.7];
% barh(Un_Model_evidence_MF_only(:,1) - Un_Model_evidence_MB_only(:,1), 'facecolor', col, 'edgecolor', 'none');
% ylim([0 size(T_Model_evidence_hybrid,1)+1])
% title('model evidence difference_dishonesty', 'FontSize', 16)
% ylabel('subject', 'FontSize', 12)
% xlabel('\Delta lme', 'FontSize', 12)
% text(-80,20,'MB', 'FontSize', 14)
% text(10,20,'MF', 'FontSize', 14)
% box off

%%
% GBF = zeros(4, 1);
% BF = zeros(30, 4);
% 
% for i = 1 : 30
%     BF(i, 1) = exp(T_Model_evidence_hybrid(i,1) - T_Model_evidence_MF_only(i,1));
%     BF(i, 2) = exp(T_Model_evidence_hybrid(i,1) - T_Model_evidence_MB_only(i,1));
%     BF(i, 3) = exp(Un_Model_evidence_hybrid(i,1) - Un_Model_evidence_MF_only(i,1));
%     BF(i, 4) = exp(Un_Model_evidence_hybrid(i,1) - Un_Model_evidence_MB_only(i,1));
% end
% 
% GBF(1,1) = prod(BF(:,1))
% GBF(2,1) = prod(BF(:,2))
% GBF(3,1) = prod(BF(:,3))
% GBF(4,1) = prod(BF(:,4))

save 'GBF.mat'