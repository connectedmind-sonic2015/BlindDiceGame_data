clear all;
tic;
%% load data 
file_list = ['subject01'; 'subject02'; 'subject03'; 'subject04'; 'subject05'; 'subject06'; 'subject07'; 'subject08'; 'subject09'; 'subject10';...
     'subject11'; 'subject12'; 'subject13'; 'subject14'; 'subject15'; 'subject16'; 'subject17'; 'subject18'; 'subject19'; 'subject20';...
     'subject21'; 'subject22'; 'subject23'; 'subject24'; 'subject25'; 'subject26'; 'subject27'; 'subject28'; 'subject29'; 'subject30']; 
 
sub_no = 30; 
ntrial = 64;
interval_size = 8;
trial = (1:1:64);

Condreward = [ 25, 50, 100, 300 ]./300;
Condpenalty = [ -25, -50, -100, -300 ]./300;
 
% Condreward = [ 25, 50, 100, 300 ];
% Condpenalty = [ -25, -50, -100, -300 ];

order = zeros(sub_no, 1); % 1 = Honesty first, 2 = Dishonesty first
T_data_summary = zeros(3, (ntrial/interval_size) ,sub_no); % 1 = openratio, 2 = mean experienced penalty, 3 = mean experienced reward
Un_data_summary = zeros(3, (ntrial/interval_size) ,sub_no); % 1 = openratio, 2 = mean experienced penalty, 3 = mean experienced reward

for sub_no = 1 : length(file_list)
    L = ls(['*' file_list(sub_no,:) '*']);
       
    for i = 1:2
        load(L(i, : )); 
        order(sub_no, 1) = B_char(1,1);
        if B_char(:, i) == 1
            Trustworthy = trial;
        else 
            Untrustworthy = trial;
        end
    end
    
    T_open = zeros( ntrial, 1 );
    T_r =zeros( ntrial, 1 );
    T_True = zeros( ntrial, 1 );
    T_trial_reward = zeros( ntrial, 1 );
    T_trial_penalty = zeros( ntrial, 1 ); 

    Un_open = zeros( ntrial, 1 );
    Un_r =zeros( ntrial, 1 );
    Un_True = zeros( ntrial, 1 );
    Un_trial_reward = zeros( ntrial, 1 );
    Un_trial_penalty = zeros( ntrial, 1 ); 

    
    for i = 1 : ntrial
        
        T_True( i, 1 ) = -(Trustworthy(i).intend - 1);
        T_trial_reward( i, 1 ) = Condreward(Trustworthy(i).uncover_reward);
        T_trial_penalty( i, 1 ) = Condpenalty(Trustworthy(i).uncover_penalty);
        T_open( i, 1 ) = Trustworthy(i).resp_U + 1;
        if T_open( i, 1 ) == 1 % 확인
            if T_True( i, 1 ) == 0 % 딜러 거짓말        
                T_r( i, 1 ) = Condreward(Trustworthy(i).uncover_reward);
            else
                T_r( i, 1 ) = Condpenalty(Trustworthy(i).uncover_penalty);
            end
        end

        Un_True( i, 1 ) = -(Untrustworthy(i).intend - 1);
        Un_trial_reward( i, 1 ) = Condreward(Untrustworthy(i).uncover_reward);
        Un_trial_penalty( i, 1 ) = Condpenalty(Untrustworthy(i).uncover_penalty);
        Un_open( i, 1 ) = Untrustworthy(i).resp_U + 1;
        if Un_open( i, 1 ) == 1 % 학인
            if Un_True( i, 1 ) == 0 % 딜러 거짓말        
                Un_r( i, 1 ) = Condreward(Untrustworthy(i).uncover_reward);
            else
                Un_r( i, 1 ) = Condpenalty(Untrustworthy(i).uncover_penalty);
            end
        end
    end
    
    
    for j = 1 : ntrial/interval_size
        T_data_summary(1, j, sub_no) = (interval_size - sum(T_open((1+interval_size*(j-1):interval_size*j),1) - 1))/interval_size;
        T_data_summary(2, j, sub_no) = mean(T_trial_penalty(find(T_open((1+interval_size*(j-1):interval_size*j),1)==1)));
        T_data_summary(3, j, sub_no) = mean(T_trial_reward(find(T_open((1+interval_size*(j-1):interval_size*j),1)==1)));
        Un_data_summary(1, j, sub_no) = (interval_size - sum(Un_open((1+interval_size*(j-1):interval_size*j),1) - 1))/interval_size;
        Un_data_summary(2, j, sub_no) = mean(Un_trial_penalty(find(Un_open((1+interval_size*(j-1):interval_size*j),1)==1)));
        Un_data_summary(3, j, sub_no) = mean(Un_trial_reward(find(Un_open((1+interval_size*(j-1):interval_size*j),1)==1)));
    end

end


