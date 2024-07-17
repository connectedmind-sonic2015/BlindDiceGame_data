% clear all

% load 'exp_GBF.mat'
F_model_evidence(:,1) = Model_evidence_hybrid;
F_model_evidence(:,2) = Model_evidence_twoalpha;
% F_model_evidence(:,2) = Model_evidence_MF_only;
% F_model_evidence(:,3) = Model_evidence_MB_only;
% P_model_evidence(:,1) = P_Model_evidence_hybrid;
% P_model_evidence(:,2) = P_Model_evidence_MF_only;
% P_model_evidence(:,3) = P_Model_evidence_MB_only;


F_alpha = ones(2,1);
num_m = 2;
num_subject = 33;

%% honesty dealer condition 
Unk = zeros(33,2);
Un = zeros(33,1);
Gnk = zeros(33,2);
beta = zeros(2,1);
cnt = 1;
diff_ratio = zeros(3,1);

while 1
    cnt = cnt + 1
    % compute gnk
    for i = 1 : num_subject
        for j = 1 : num_m
            Unk(i, j) = exp (F_model_evidence(i, j) + psi(F_alpha(j, (cnt - 1))) - psi(sum(F_alpha(:, (cnt - 1 )))));
        end
        Un(i, 1) = sum(Unk(i,:));
    end
    
    % compute β
    for i = 1 : num_subject
        for j = 1 : num_m
            Gnk(i, j) = Unk(i, j)/Un(i, 1);
        end
    end

    for j = 1 : num_m
        beta(j, 1) = sum(Gnk(:, j));
    end

    % update α = α0 + β
    F_alpha( : , cnt ) = F_alpha( :, cnt-1) + beta;
    diff_ratio = F_alpha( :, cnt - 1) ./ F_alpha( :, cnt )
    if diff_ratio >= 0.99  % until convergence
%     diff = T_alpha( :, cnt - 1) - T_alpha( :, cnt );
%     if abs(diff) <= 0.00001  % until convergence
       F_alpha( :, cnt)
        break
    end
end

% compute Rk,exceedance probability


%% dishonesty dealer condition
% Unk = zeros(33,3);
% Un = zeros(33,1);
% Gnk = zeros(33,3);
% beta = zeros(3,1);
% cnt = 1;
% diff_ratio = zeros(3,1);
% 
% while 1
%     cnt = cnt + 1
%     % compute gnk
%     for i = 1 : num_subject
%         for j = 1 : num_m
%             Unk(i, j) = exp (P_model_evidence(i, j) + psi(P_alpha(j, (cnt - 1))) - psi(sum(P_alpha(:, (cnt - 1)))));
%         end
%         Un(i, 1) = sum(Unk(i,:));
%     end
%     
%     % compute β
%     for i = 1 : num_subject
%         for j = 1 : num_m
%             Gnk(i, j) = Unk(i, j)/Un(i, 1);
%         end
%     end
% 
%     for j = 1 : num_m
%         beta(j, 1) = sum(Gnk(:, j));
%     end
% 
%     % update α = α0 + β
%     P_alpha( : , cnt ) = P_alpha( :, cnt-1) + beta;
%     diff_ratio = P_alpha( :, cnt - 1) ./ P_alpha( :, cnt )
%     if diff_ratio >= 0.99  % until convergence
% %     diff = UN_alpha( :, cnt - 1) - UN_alpha( :, cnt )
% %     if abs(diff) <= 0.00001  % until convergence
%        P_alpha( :, cnt)
%         break
%     end
% end