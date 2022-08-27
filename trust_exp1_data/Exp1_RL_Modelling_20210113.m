clear all;
tic;
%% load data 
file_list = ['subject01'; 'subject02'; 'subject03'; 'subject04'; 'subject05'; 'subject06'; 'subject07'; 'subject08'; 'subject09'; 'subject10';...
     'subject11'; 'subject12'; 'subject13'; 'subject14'; 'subject15'; 'subject16'; 'subject17'; 'subject18'; 'subject19'; 'subject20';...
     'subject21'; 'subject22'; 'subject23'; 'subject24'; 'subject25'; 'subject26'; 'subject27'; 'subject28'; 'subject29'; 'subject30']; 
%

sub_no = 29; 
% ntrial = 64;
ntrial = 20;
trial = (1:1:64);

Condreward = [ 25, 50, 100, 300 ]./300;
Condpenalty = [ -25, -50, -100, -300 ]./300;


beta_lb = 0; 
rho_lb = 0;
beta_ub = 10; 
rho_ub = 10;
options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set', 'FunValCheck', 'on');
MAP =zeros(1,1);

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
            % model 2alphas
            % Trustworthy
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
%             x0 = [rand, exprnd(1), rand, exprnd(1), exprnd(1), rand];
            x0 = [rand, rand, exprnd(1), rand];
            dealer = 1; %신뢰

            [Xfit, T_NegMAP] = fmincon(@(x) computeMAP_2alpha (x, T_a, T_r, T_True, T_trial_reward, T_trial_penalty), x0, [], [], [], [],...
                [0, 0, beta_lb, 0], [1, 1, beta_ub, 1], [], options);
            T_pars_MAP(sub_no, :) = [Xfit, T_NegMAP];
            T_MAP(sub_no, 1) = T_NegMAP;
            
            % Untrustworthy
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
%             x0 = [rand, exprnd(1), rand, exprnd(1), exprnd(1), rand];
            x0 = [rand, rand, exprnd(1), rand];
            dealer = 0; %비신뢰

            [Xfit, UN_NegMAP] = fmincon(@(x) computeMAP_2alpha (x, Un_a, Un_r, Un_True, Un_trial_reward, Un_trial_penalty), x0, [], [], [], [],...
                [0, 0, beta_lb, 0], [1, 1, beta_ub, 1], [], options);
            Un_pars_MAP(sub_no, :) = [Xfit, UN_NegMAP];
            UN_MAP(sub_no, 1) = UN_NegMAP;


    end


%% compute LLH

% Model_5

function NegMAP = computeMAP_2alpha (p1, a, r,  True,  trial_reward, trial_penalty)

    alpha_MF = p1(1);
    alpha_MB = p1(2);
    beta = p1(3);
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

