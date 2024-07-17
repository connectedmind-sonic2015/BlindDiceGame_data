clear all;
tic;

%% AIC: −2 × log-likelihood + 2 × number of parameters
%% BIC: −2 × log-likelihood + log(Number of observations) × number of parameters

%% load data 
%  'subject01' 45trials; 'subject02'45trials; 'subject17'돈떨어져 P 44; 'subject32' 돈떨어져 P 43; 'subject34' 파일 비어있음.; 'subject54' 돈떨어져 P 42;
file_list = ['subject01'; 'subject02'; 'subject03'; 'subject04'; 'subject05'; 'subject06'; 'subject07'; 'subject08'; 'subject09'; 'subject10';...
     'subject11'; 'subject12'; 'subject13'; 'subject14'; 'subject15'; 'subject16'; 'subject17'; 'subject18'; 'subject19'; 'subject20';...
     'subject21'; 'subject22'; 'subject23'; 'subject24'; 'subject25'; 'subject26'; 'subject27'; 'subject28'; 'subject29'; 'subject30';...
     'subject31'; 'subject32'; 'subject33'];


sub_num = length(file_list);
beta_lb = 0; 
rho_lb = 0;
beta_ub = 10; 
rho_ub = 10;
options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set', 'FunValCheck', 'on');


ex2_Hybrid_pars_MAP = zeros(sub_num, 1, 5);


%% AIC: −2 × log-likelihood + 2 × number of parameters
%% BIC: −2 × log-likelihood + log(Number of observations) × number of parameters


cur_dir = pwd;
for sub_num = 1 : length(file_list)
    cd (cur_dir)

    target_dir = file_list(sub_num,:);
    cd (target_dir)

%

    Monetary_data = deblank(ls(['*BDG_financial_2020*']));
    Mon_data = readtable(Monetary_data);
    M_data = table2struct(Mon_data);
    M_ntrial = length(M_data);
    M_trial = (1:1:M_ntrial);


% load variables from monetary/peg data
    M_honesty =zeros( M_ntrial, 1 );
    M_open_reward = zeros( M_ntrial, 1 );
    M_open_penalty = zeros( M_ntrial, 1 );
    M_choice = zeros( M_ntrial, 1 );
    M_choice_reward = zeros( M_ntrial, 1 ); 

%   
    M_count = 0;


    for i = 1 : M_ntrial
        if string(M_data(i).whether_to_lie) == "truth"
            M_honesty(i, 1) = 1; %진실
        elseif string(M_data(i).whether_to_lie) == "lie"
            M_honesty(i, 1) = 0; %거짓
        end
        M_trial_rating( i, 1 ) = M_data(i).trust_rating;
        M_open_reward(i, 1) = (M_data(i).open_reward/1000)*3; % 1000
        M_open_penalty(i, 1) = (M_data(i).open_penalty/1000)*3; % 1000
        if string(M_data(i).user_response_open) == 'open'
            M_choice(i, 1) = 1; % 확인
            M_choice_reward(i, 1) = (M_data(i).choice_reward/1000)*3; % 1000
            if string(M_data(i).whether_to_lie) == "truth"
                M_count = M_count + 1;
            end
        elseif string(M_data(i).user_response_open) == 'not_open'
            M_choice(i, 1) = 2; % 확인 안함
            M_choice_reward(i, 1) = 0;
        end
    end


%% parameter fitting
    for j = 1: M_ntrial
    % hybrid model with 2 alphas
    % fianancial
        options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
        x0 = [rand, rand, rand, exprnd(1)]; % MF alpha, MB alpha, w, t0, beta
    
        [Xfit, T_NegMAP,exitflag,output,lambda,grad,hessian] = fmincon(@(x) computeMAP (x, M_choice(1:j, 1), M_choice_reward(1:j, 1), M_honesty(1:j, 1), M_open_reward(1:j, 1), M_open_penalty(1:j, 1)), x0, [], [], [], [],...
            [0, 0, 0, beta_lb,], [1, 1, 1, beta_ub], [], options);
        ex2_Hybrid_pars_MAP(sub_num, j, :) = [Xfit, T_NegMAP];
    
    end

end

 


%% compute MAP


function NegMAP = computeMAP (p1, a, r,  True,  trial_reward, trial_penalty)

    alpha_MF = p1(1);
    alpha_MB = p1(2);
    w = p1(3);
    t0 = 0.5;
    beta = p1(4);

    T = length(a);

    QMF = [0, 0];
    netQ = [0, 0];
    deltaMF = 0;
    deltaMB = 0;
    TMB = zeros(T, 1);
    TMB(1,1) = t0;
    QMB = [0, 0]; 
    
    % loop over all trial
    for t = 1:T
        QMB(1, 1)= ((1 - TMB(t)) * trial_reward(t)) + (TMB(t) * trial_penalty(t));
      
        netQ = ((1 - w) * QMF) + (w * QMB);
        p = exp(beta*netQ) / sum(exp(beta*netQ));

        if p(a(t)) == 0
            choiceProb(t) = 1;% 1 이 아닌 아주 작은 수??
        else
            choiceProb(t) = p(a(t));
        end
        
        % update values
        if a(t) == 1
            
            deltaMF = r(t) - QMF(a(t));
            QMF(a(t)) = QMF(a(t)) + alpha_MF * deltaMF;

            deltaMB = True(t) - TMB(t);
            if t < T
                TMB(t + 1) = TMB(t) + alpha_MB * deltaMB;
            end
            TMB(t);         
        else
            if t < T
               TMB(t + 1) = TMB(t);
            end
         end
    end

    % compute negative log-likelihood
   % NegLL = -sum(log(choiceProb));
%      NegMAP = -sum(log(choiceProb)); 
%      
    alpha_MFProb = betapdf(alpha_MF, 1.1, 1.1);
    alpha_MBProb = betapdf(alpha_MB, 1.1, 1.1);    
    wProb = betapdf(w, 1.1, 1.1);
    betaProb = gampdf(beta, 1.4, 5);
    
    NegMAP = -sum(log(choiceProb) + log(alpha_MFProb) + log(alpha_MBProb) + log(wProb) + log(betaProb) );
 
end
