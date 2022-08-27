clear all;
tic;
%% load data 
file_list = ['subject01'; 'subject02'; 'subject03'; 'subject04'; 'subject05'; 'subject06'; 'subject07'; 'subject08'; 'subject09'; 'subject10';...
     'subject11'; 'subject12'; 'subject13'; 'subject14'; 'subject15'; 'subject16'; 'subject17'; 'subject18'; 'subject19'; 'subject20';...
     'subject21'; 'subject22'; 'subject23'; 'subject24'; 'subject25'; 'subject26'; 'subject27'; 'subject28'; 'subject29'; 'subject30']; 
 
sub_no = 30; 
ntrial = 64;
trial = (1:1:64);

Condreward = [ 25, 50, 100, 300 ]./300;
Condpenalty = [ -25, -50, -100, -300 ]./300;

beta_lb = 0; 
beta_ub = 10; 

options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set', 'FunValCheck', 'on');

Hybrid_T_MAP = zeros(1,1);
Hybrid_UN_MAP =zeros(1,1);
Hybrid_T_hessian = zeros(4,4,sub_no);
Hybrid_UN_hessian = zeros(4,4,sub_no);

MFonly_T_MAP = zeros(1,1);
MFonly_UN_MAP =zeros(1,1);
MFonly_T_hessian = zeros(3,3,sub_no);
MFonly_UN_hessian = zeros(3,3,sub_no);

MBonly_T_MAP = zeros(1,1);
MBonly_UN_MAP =zeros(1,1);
MBonly_T_hessian = zeros(3,3,sub_no);
MBonly_UN_hessian = zeros(3,3,sub_no);


betaCons_T_MAP = zeros(1,1);
betaCons_UN_MAP =zeros(1,1);
betaCons_T_hessian = zeros(3,3,sub_no);
betaCons_UN_hessian = zeros(3,3,sub_no);

T_Model_evidence_hybrid = zeros(sub_no, 1);
T_Model_evidence_MF_only = zeros(1,sub_no, 1);
T_Model_evidence_MB_only = zeros(1,sub_no, 1);

Un_Model_evidence_hybrid = zeros(sub_no, 1);
Un_Model_evidence_MF_only = zeros(sub_no, 1);
Un_Model_evidence_MB_only = zeros(sub_no, 1);

for sub_no = 1 : length(file_list)
    L = ls(['*' file_list(sub_no,:) '*']);
   
    
    for i = 1:2
        load(L(i, : ));
        if B_char(:, i) == 1
            Trustworthy = trial;
        else 
            Untrustworthy = trial;
        end
    end
    
    T_a = zeros( ntrial, 1 );
    T_r =zeros( ntrial, 1 );
    T_True = zeros( ntrial, 1 );
    T_trial_reward = zeros( ntrial, 1 );
    T_trial_penalty = zeros( ntrial, 1 ); 

    Un_a = zeros( ntrial, 1 );
    Un_r =zeros( ntrial, 1 );
    Un_True = zeros( ntrial, 1 );
    Un_trial_reward = zeros( ntrial, 1 );
    Un_trial_penalty = zeros( ntrial, 1 ); 

    
    for i = 1 : ntrial
        
        T_True( i, 1 ) = -(Trustworthy(i).intend - 1);
        T_trial_reward( i, 1 ) = Condreward(Trustworthy(i).uncover_reward);
        T_trial_penalty( i, 1 ) = Condpenalty(Trustworthy(i).uncover_penalty);
        T_a( i, 1 ) = Trustworthy(i).resp_U + 1;
        if T_a( i, 1 ) == 1 % 확인
            if T_True( i, 1 ) == 0 % 딜러 거짓말        
                T_r( i, 1 ) = Condreward(Trustworthy(i).uncover_reward);
            else
                T_r( i, 1 ) = Condpenalty(Trustworthy(i).uncover_penalty);
            end
        end

        Un_True( i, 1 ) = -(Untrustworthy(i).intend - 1);
        Un_trial_reward( i, 1 ) = Condreward(Untrustworthy(i).uncover_reward);
        Un_trial_penalty( i, 1 ) = Condpenalty(Untrustworthy(i).uncover_penalty);
        Un_a( i, 1 ) = Untrustworthy(i).resp_U + 1;
        if Un_a( i, 1 ) == 1 % 학인
            if Un_True( i, 1 ) == 0 % 딜러 거짓말        
                Un_r( i, 1 ) = Condreward(Untrustworthy(i).uncover_reward);
            else
                Un_r( i, 1 ) = Condpenalty(Untrustworthy(i).uncover_penalty);
            end
        end
    end

        
% initial setting

% parameter fitting
            % hybrid model 2alphas
            % Trustworthy
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
            x0 = [rand, rand, exprnd(1), rand];
            dealer = 1; %신뢰

            [Xfit, T_NegMAP,exitflag,output,lambda,grad,hessian] = fmincon(@(x) computeMAP (x, T_a, T_r, T_True, T_trial_reward, T_trial_penalty), x0, [], [], [], [],...
                [0, 0, beta_lb, 0], [1, 1, beta_ub, 1], [], options);
            Hybrid_T_pars_MAP(sub_no, :) = [Xfit, T_NegMAP];
            Hybrid_T_MAP(sub_no, 1) = T_NegMAP;
            Hybrid_T_hessian(: ,: ,sub_no) = hessian;
            
            % hybrid Untrustworthy
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
            x0 = [rand, rand, exprnd(1), rand];
            dealer = 0; %비신뢰

            [Xfit, UN_NegMAP, exitflag,output,lambda,grad,hessian] = fmincon(@(x) computeMAP (x, Un_a, Un_r, Un_True, Un_trial_reward, Un_trial_penalty), x0, [], [], [], [],...
                [0, 0, beta_lb, 0], [1, 1, beta_ub, 1], [], options);
            Hybrid_Un_pars_MAP(sub_no, :) = [Xfit, UN_NegMAP];
            Hybrid_UN_MAP(sub_no, 1) = UN_NegMAP;
            Hybrid_UN_hessian(:, :, sub_no) = hessian;

            %% MF only
            % Trustworthy
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
            x0 = [rand, rand, exprnd(1)];
            dealer = 1; %신뢰

            [Xfit, T_NegMAP,exitflag,output,lambda,grad,hessian] = fmincon(@(x) computeMAP_MFonly (x, T_a, T_r, T_True, T_trial_reward, T_trial_penalty), x0, [], [], [], [],...
                [0, 0, beta_lb], [1, 1, beta_ub], [], options);
            MFonly_T_pars_MAP(sub_no, :) = [Xfit, T_NegMAP];
            MFonly_T_MAP(sub_no, 1) = T_NegMAP;
            MFonly_T_hessian(: ,: ,sub_no) = hessian;
            
            % hybrid Untrustworthy
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
            x0 = [rand, rand, exprnd(1)];
            dealer = 0; %비신뢰

            [Xfit, UN_NegMAP, exitflag,output,lambda,grad,hessian] = fmincon(@(x) computeMAP_MFonly (x, Un_a, Un_r, Un_True, Un_trial_reward, Un_trial_penalty), x0, [], [], [], [],...
                [0, 0, beta_lb], [1, 1, beta_ub], [], options);
            MFonly_Un_pars_MAP(sub_no, :) = [Xfit, UN_NegMAP];
            MFonly_UN_MAP(sub_no, 1) = UN_NegMAP;
            MFonly_UN_hessian(:, :, sub_no) = hessian;

            %% Model Based only
            % Trustworthy
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
            x0 = [rand, rand, exprnd(1)];
            dealer = 1; %신뢰

            [Xfit, T_NegMAP,exitflag,output,lambda,grad,hessian] = fmincon(@(x) computeMAP_MBonly (x, T_a, T_r, T_True, T_trial_reward, T_trial_penalty), x0, [], [], [], [],...
                [0, 0, beta_lb], [1, 1, beta_ub], [], options);
            MBonly_T_pars_MAP(sub_no, :) = [Xfit, T_NegMAP];
            MBonly_T_MAP(sub_no, 1) = T_NegMAP;
            MBonly_T_hessian(: ,: ,sub_no) = hessian;
            
            % Untrustworthy
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
            x0 = [rand, rand, exprnd(1)];
            dealer = 0; %비신뢰

            [Xfit, UN_NegMAP, exitflag,output,lambda,grad,hessian] = fmincon(@(x) computeMAP_MBonly (x, Un_a, Un_r, Un_True, Un_trial_reward, Un_trial_penalty), x0, [], [], [], [],...
                [0, 0, beta_lb], [1, 1, beta_ub], [], options);
            MBonly_Un_pars_MAP(sub_no, :) = [Xfit, UN_NegMAP];
            MBonly_UN_MAP(sub_no, 1) = UN_NegMAP;
            MBonly_UN_hessian(:, :, sub_no) = hessian;
 %% beta = 1 
            % Trustworthy
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
            x0 = [rand, rand, rand];
            dealer = 1; %신뢰

            [Xfit, T_NegMAP,exitflag,output,lambda,grad,hessian] = fmincon(@(x) computeMAP_betaCons (x, T_a, T_r, T_True, T_trial_reward, T_trial_penalty), x0, [], [], [], [],...
                [0, 0, 0], [1, 1, 1], [], options);
            betaCons_T_pars_MAP(sub_no, :) = [Xfit, T_NegMAP];
            betaCons_T_MAP(sub_no, 1) = T_NegMAP;
            betaCons_T_hessian(: ,: ,sub_no) = hessian;
            
            % Untrustworthy
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
            x0 = [rand, rand, rand];
            dealer = 0; %비신뢰

            [Xfit, UN_NegMAP, exitflag,output,lambda,grad,hessian] = fmincon(@(x) computeMAP_betaCons (x, Un_a, Un_r, Un_True, Un_trial_reward, Un_trial_penalty), x0, [], [], [], [],...
                [0, 0, 0], [1, 1, 1], [], options);
            betaCons_Un_pars_MAP(sub_no, :) = [Xfit, UN_NegMAP];
            betaCons_UN_MAP(sub_no, 1) = UN_NegMAP;
            betaCons_UN_hessian(:, :, sub_no) = hessian;
            
end

for k = 1 : sub_no
    T_Model_evidence_hybrid(k, 1) = Hybrid_T_MAP(k, 1) + 4/2*log(pi) + 1/2*log(det(Hybrid_T_hessian(:, :, k)));
    T_Model_evidence_MF_only(k, 1) = MFonly_T_MAP(k, 1) + 3/2*log(pi) + 1/2*log(det(MFonly_T_hessian(:, :, k)));
    T_Model_evidence_MB_only(k, 1) = MBonly_T_MAP(k, 1) + 3/2*log(pi) + 1/2*log(det(MBonly_T_hessian(:, :, k)));    

    Un_Model_evidence_hybrid(k, 1) = Hybrid_UN_MAP(k, 1) + 4/2*log(pi) + 1/2*log(det(Hybrid_UN_hessian(:, :, k)));
    Un_Model_evidence_MF_only(k, 1) = MFonly_UN_MAP(k, 1) + 3/2*log(pi) + 1/2*log(det(MFonly_UN_hessian(:, :, k)));
    Un_Model_evidence_MB_only(k, 1) = MBonly_UN_MAP(k, 1) + 3/2*log(pi) + 1/2*log(det(MBonly_UN_hessian(:, :, k)));    
end



%% compute LLH


function NegMAP = computeMAP (p1, a, r,  True,  trial_reward, trial_penalty)

    alpha_MF = p1(1);
    alpha_MB = p1(2);
    beta = p1(3);
 %   p1(4) =1;
    w = p1(4);

    T = length(a);

    QMF = [0, 0];
    netQ = [0, 0];
    deltaMF = 0;
    deltaMB = 0;
    TMB = zeros(T, 1);
    TMB(1,1) = 0.5;
    QMB = 0; 
    
    % loop over all trial
    for t = 1:T
        QMB = ((1 - TMB(t)) * trial_reward(t)) + (TMB(t) * trial_penalty(t));
        netQ(a(t)) = ((1 - w) * QMF(a(t))) + (w * QMB);
        p = exp(beta*netQ) / sum(exp(beta*netQ));

        if p(a(t)) == 0
            choiceProb(t) = 1;% 1 이 아닌 아주 작은 수??
        else
            choiceProb(t) = p(a(t));
        end
        
        % update values
        if a(t) == 1
            
            deltaMF = r(t) - QMF(a(t));
            QMF(a(t)) = QMF(a(t)) + alpha_MF * deltaMF

            deltaMB = True(t) - TMB(t);
            if t < T
                TMB(t + 1) = TMB(t) + alpha_MB * deltaMB;
            end
            TMB(t)           
        else
            if t < T
               TMB(t + 1) = TMB(t);
            end
         end
    end

    % compute negative log-likelihood
    NegLL = -sum(log(choiceProb));
        
    alpha_MFProb = betapdf(alpha_MF, 1.1, 1.1);
    alpha_MBProb = betapdf(alpha_MB, 1.1, 1.1);
    betaProb = gampdf(beta, 1.2, 5);
    wProb = betapdf(w, 1.1, 1.1);    
    
    NegMAP = -sum(log(choiceProb) + log(alpha_MFProb) + log(alpha_MBProb)  + log(betaProb) + log(wProb));
 
end

function NegMAP = computeMAP_MFonly (p1, a, r,  True,  trial_reward, trial_penalty)

    alpha_MF = p1(1);
    alpha_MB = p1(2);
    beta = p1(3);
    w = 0;

    T = length(a);

    QMF = [0, 0];
    netQ = [0, 0];
    deltaMF = 0;
    deltaMB = 0;
    TMB = zeros(T, 1);
    TMB(1,1) = 0.5;
    QMB = 0; 
    
    % loop over all trial
    for t = 1:T
        QMB = ((1 - TMB(t)) * trial_reward(t)) + (TMB(t) * trial_penalty(t));
        netQ(a(t)) = ((1 - w) * QMF(a(t))) + (w * QMB);
        p = exp(beta*netQ) / sum(exp(beta*netQ));

        if p(a(t)) == 0
            choiceProb(t) = 1;% 1 이 아닌 아주 작은 수??
        else
            choiceProb(t) = p(a(t));
        end
        
        % update values
        if a(t) == 1
            
            deltaMF = r(t) - QMF(a(t));
            QMF(a(t)) = QMF(a(t)) + alpha_MF * deltaMF

            deltaMB = True(t) - TMB(t);
            if t < T
                TMB(t + 1) = TMB(t) + alpha_MB * deltaMB;
            end
            TMB(t)           
        else
            if t < T
               TMB(t + 1) = TMB(t);
            end
         end
    end

    % compute negative log-likelihood
    NegLL = -sum(log(choiceProb));
        
    alpha_MFProb = betapdf(alpha_MF, 1.1, 1.1);
    alpha_MBProb = betapdf(alpha_MB, 1.1, 1.1);
    betaProb = gampdf(beta, 1.2, 5);
%     wProb = betapdf(w, 1.1, 1.1);    
    
%     NegMAP = -sum(log(choiceProb) + log(alpha_MFProb) + log(alpha_MBProb)  + log(betaProb) + log(wProb));
    NegMAP = -sum(log(choiceProb) + log(alpha_MFProb) + log(alpha_MBProb)  + log(betaProb));
 
end

function NegMAP = computeMAP_MBonly (p1, a, r,  True,  trial_reward, trial_penalty)

    alpha_MF = p1(1);
    alpha_MB = p1(2);
    beta = p1(3);
    w = 1;

    T = length(a);

    QMF = [0, 0];
    netQ = [0, 0];
    deltaMF = 0;
    deltaMB = 0;
    TMB = zeros(T, 1);
    TMB(1,1) = 0.5;
    QMB = 0; 
    
    % loop over all trial
    for t = 1:T
        QMB = ((1 - TMB(t)) * trial_reward(t)) + (TMB(t) * trial_penalty(t));
        netQ(a(t)) = ((1 - w) * QMF(a(t))) + (w * QMB);
        p = exp(beta*netQ) / sum(exp(beta*netQ));

        if p(a(t)) == 0
            choiceProb(t) = 1;% 1 이 아닌 아주 작은 수??
        else
            choiceProb(t) = p(a(t));
        end
        
        % update values
        if a(t) == 1
            
            deltaMF = r(t) - QMF(a(t));
            QMF(a(t)) = QMF(a(t)) + alpha_MF * deltaMF

            deltaMB = True(t) - TMB(t);
            if t < T
                TMB(t + 1) = TMB(t) + alpha_MB * deltaMB;
            end
            TMB(t)           
        else
            if t < T
               TMB(t + 1) = TMB(t);
            end
         end
    end

    % compute negative log-likelihood
    NegLL = -sum(log(choiceProb));
        
    alpha_MFProb = betapdf(alpha_MF, 1.1, 1.1);
    alpha_MBProb = betapdf(alpha_MB, 1.1, 1.1);
    betaProb = gampdf(beta, 1.2, 5);
%     wProb = betapdf(w, 1.1, 1.1);    
    
%     NegMAP = -sum(log(choiceProb) + log(alpha_MFProb) + log(alpha_MBProb)  + log(betaProb) + log(wProb));
    NegMAP = -sum(log(choiceProb) + log(alpha_MFProb) + log(alpha_MBProb)  + log(betaProb));
 
end

function NegMAP = computeMAP_betaCons (p1, a, r,  True,  trial_reward, trial_penalty)

    alpha_MF = p1(1);
    alpha_MB = p1(2);
    w = p1(3);
    beta = 1;

    T = length(a);

    QMF = [0, 0];
    netQ = [0, 0];
    deltaMF = 0;
    deltaMB = 0;
    TMB = zeros(T, 1);
    TMB(1,1) = 0.5;
    QMB = 0; 
    
    % loop over all trial
    for t = 1:T
        QMB = ((1 - TMB(t)) * trial_reward(t)) + (TMB(t) * trial_penalty(t));
        netQ(a(t)) = ((1 - w) * QMF(a(t))) + (w * QMB);
        p = exp(beta*netQ) / sum(exp(beta*netQ));

        if p(a(t)) == 0
            choiceProb(t) = 1;% 1 이 아닌 아주 작은 수??
        else
            choiceProb(t) = p(a(t));
        end
        
        % update values
        if a(t) == 1
            
            deltaMF = r(t) - QMF(a(t));
            QMF(a(t)) = QMF(a(t)) + alpha_MF * deltaMF

            deltaMB = True(t) - TMB(t);
            if t < T
                TMB(t + 1) = TMB(t) + alpha_MB * deltaMB;
            end
            TMB(t)           
        else
            if t < T
               TMB(t + 1) = TMB(t);
            end
         end
    end

    % compute negative log-likelihood
    NegLL = -sum(log(choiceProb));
        
    alpha_MFProb = betapdf(alpha_MF, 1.1, 1.1);
    alpha_MBProb = betapdf(alpha_MB, 1.1, 1.1);
%     betaProb = gampdf(beta, 1.2, 5);
    wProb = betapdf(w, 1.1, 1.1);    
    
    NegMAP = -sum(log(choiceProb) + log(alpha_MFProb) + log(alpha_MBProb) + log(wProb));
 
end

