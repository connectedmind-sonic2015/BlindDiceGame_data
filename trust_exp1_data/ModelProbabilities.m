clear all

load 'modelfiiting_evidence.mat'
T_model_evidence(:,1) = T_Model_evidence_hybrid;
T_model_evidence(:,2) = T_Model_evidence_MF_only;
T_model_evidence(:,3) = T_Model_evidence_MB_only;
UN_model_evidence(:,1) = Un_Model_evidence_hybrid;
UN_model_evidence(:,2) = Un_Model_evidence_MF_only;
UN_model_evidence(:,3) = Un_Model_evidence_MB_only;

T_alpha = ones(3,1);
UN_alpha = ones(3,1);
num_m = 3;
num_subject = 30;

%% honesty dealer condition 
Unk = zeros(30,3);
Un = zeros(30,1);
Gnk = zeros(30,3);
beta = zeros(3,1);
cnt = 1;

while 1
    cnt = cnt + 1
    % compute gnk
    for i = 1 : num_subject
        for j = 1 : num_m
            Unk(i, j) = exp (T_model_evidence(i, j) + psi(T_alpha(j, 1)) - psi(sum(T_alpha(:, 1))));
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
    T_alpha( : , cnt ) = T_alpha( :, cnt-1) + beta;
    diff_ratio = T_alpha( :, cnt - 1) ./ T_alpha( :, cnt )
    if mean(diff_ratio) >= 0.99999  % until convergence
       T_alpha( :, cnt)
        break
    end
end

% compute Rk,exceedance probability


%% dishonesty dealer condition
Unk = zeros(30,3);
Un = zeros(30,1);
Gnk = zeros(30,3);
beta = zeros(3,1);
cnt = 1;

while 1
    cnt = cnt + 1
    % compute gnk
    for i = 1 : num_subject
        for j = 1 : num_m
            Unk(i, j) = exp (UN_model_evidence(i, j) + psi(UN_alpha(j, 1)) - psi(sum(UN_alpha(:, 1))));
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
    UN_alpha( : , cnt ) = UN_alpha( :, cnt-1) + beta;
    diff_ratio = UN_alpha( :, cnt - 1) ./ UN_alpha( :, cnt )
    if mean(diff_ratio) >= 0.99999  % until convergence
       UN_alpha( :, cnt)
        break
    end
end