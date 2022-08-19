clear all;
tic;
%% load data 
file_list = ['subject01'; 'subject02'; 'subject03'; 'subject04'; 'subject07'; 'subject08'; 'subject09'; 'subject10';...
     'subject11'; 'subject12'; 'subject13'; 'subject14'; 'subject15'; 'subject16'; 'subject17'; 'subject18'; 'subject19'; 'subject20';...
     'subject21'; 'subject22'; 'subject23'; 'subject24'; 'subject25'; 'subject26'; 'subject27'; 'subject28'; 'subject29'; 'subject30'; 'subject31'; 'subject32'; 'subject33']; 
 
sub_num = 31; 
ntrial = 36;
trial = (1:1:36);

dealer_hoenty = zeros( sub_num, 1 );
experieced_dealer_hoenty = zeros( sub_num, 1 );
% Condreward = [ 200, 500, 1000 ]./1000;
% Condpenalty = [ -200, -500, -1000 ]./1000;

% pars_MAP = zeros(sub_no, ntrial, 5);


beta_lb = 0; 
rho_lb = 0;
beta_ub = 10; 
rho_ub = 10;
options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set', 'FunValCheck', 'on');


for sub_no = 1 : length(file_list)
    L = ls(['*' file_list(sub_no,:) '*']);
    
    opts = detectImportOptions(L);
    varNames = opts.VariableNames ; % variable names
    varTypes = {'double', 'double', 'categorical','categorical','categorical','categorical','double',...   
                'double','double','double','categorical','double','double'}; 
    opts = setvartype(opts,varNames,varTypes);
%     opts = setvartype(opts,varNames,varTypes, 'PreserveVariableNames', true);
    
    T = readtable(L,opts);
%     T = readtable(L, 'PreserveVariableNames', true);
%     T = readtable('input_table.xlsx','PreserveVariableNames',true);

    block= table2struct(T);
    
    
    open = zeros( ntrial, 1 );
    reward =zeros( ntrial, 1 );
    True = zeros( ntrial, 1 );
    trial_reward = zeros( ntrial, 1 );
    trial_penalty = zeros( ntrial, 1 ); 
    trial_rating = zeros( ntrial, 1 ); 

    count = 0;
    for i = 1 : ntrial
        True( i, 1 ) = block(i).whether_to_lie == 'truth';% 딜러 거짓말 여부 1 진실, 0 거짓
        trial_rating( i, 1 ) = block(i).trust_rating;
        trial_reward( i, 1 ) = block(i).open_reward/1000;
        trial_penalty( i, 1 ) = block(i).open_penalty/1000;
        if block(i).user_response_open == 'open'; % 피험자 open 행동 여부 1 확인, 2 확인 안함  
            open( i, 1 ) = 1;
        else
            open( i, 1 ) = 2;
        end
        if open( i, 1 ) == 1 % 확인
            if True( i, 1 ) == 0 % 딜러 거짓말        
                reward( i, 1 ) = block(i).open_reward/1000;
            else
                reward( i, 1 ) = block(i).open_penalty/1000;
                count = count + 1; 
            end
        end    
    end
    
    dealer_hoenty(sub_no, 1) = sum(True)/ntrial;
    experieced_dealer_hoenty(sub_no, 1) = count /length(find(open == 1)); % 5/open한 경우 dealer가 진실을 말한 횟수
    mean_rating(sub_no, 1) = mean(trial_rating( :, 1 ));
 
%     for j = 1 : ntrial
% initial setting

% parameter fitting
            % model 2alphas
            options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
%             x0 = [rand, exprnd(1), rand, exprnd(1), exprnd(1), rand];
            x0 = [rand, rand, exprnd(1), rand];
%             dealer = 1; %신뢰

            [Xfit, NegMAP] = fmincon(@(x) computeMAP_2alpha (x, open, reward, True, trial_reward, trial_penalty), x0, [], [], [], [],...
                [0, 0, beta_lb, 0], [1, 1, beta_ub, 1], [], options);
%             pars_MAP(sub_no, j, :) = [Xfit, NegMAP];
            pars_MAP(sub_no, :) = [Xfit, NegMAP];
            

%     end
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
%         QMB = (TMB(t) * trial_penalty(t));
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
    NegLL = -sum(log(choiceProb));
        
    alpha_MFProb = betapdf(alpha_MF, 1.1, 1.1);
    alpha_MBProb = betapdf(alpha_MB, 1.1, 1.1);
    betaProb = gampdf(beta, 1.2, 5);
    wProb = betapdf(w, 1.1, 1.1);    
    
     NegMAP = -sum(log(choiceProb) + log(alpha_MFProb) + log(alpha_MBProb)  + log(betaProb) + log(wProb));
 
end

