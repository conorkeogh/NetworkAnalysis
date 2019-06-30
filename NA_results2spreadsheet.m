%% Convert results arrays into 2D spreadsheets for analysis
% Copy/paste from results2spreadsheet.m
% Adapted for responders

% i.e. to send to stats guy to do stuff in R or whatever... Might also make
% future analysis easier or something?

% Load power in each band - subject x trial x electrode x band
load('/Volumes/File Storage/JoVE/Data/powerbands.mat');

% Set current position to beginning of spreadsheet
currentline = 1;

% Cycle through subjects
for subj = 1:size(avg,1)
    
    % Cycle through each trial
    for trial = 1:size(avg,2)
        
        % Check if current trial exists
        if sum(sum(avg(subj,trial,:,:))) > 0 % All values will be 0 if doesn't exist
            
            % Set up variables to add to sheet (not necessary to assign
            % specific variables, but easier to re-organise later if
            % necessary...)
            Rfrontal = squeeze(avg(subj, trial, 1, :))';
            Rparietal = squeeze(avg(subj, trial, 2, :))';
            Roccipital = squeeze(avg(subj, trial, 3, :))';
            Rtemporal = squeeze(avg(subj, trial, 4, :))';
            Lfrontal = squeeze(avg(subj, trial, 5, :))';
            Lparietal = squeeze(avg(subj, trial, 6, :))';
            Ltemporal = squeeze(avg(subj, trial, 7, :))';
            Loccipital = squeeze(avg(subj, trial, 8, :))';
            
            overall = squeeze(mean(avg(subj, trial, :, :), 3))';
            Rhemi = squeeze(mean(avg(subj, trial, 1:4, :), 3))';
            Lhemi = squeeze(mean(avg(subj, trial, 5:8, :), 3))';
            
            % Repeat for standard deviations:
            Rfrontal_std = squeeze(std_all(subj, trial, 1, :))';
            Rparietal_std = squeeze(std_all(subj, trial, 2, :))';
            Roccipital_std = squeeze(std_all(subj, trial, 3, :))';
            Rtemporal_std = squeeze(std_all(subj, trial, 4, :))';
            Lfrontal_std = squeeze(std_all(subj, trial, 5, :))';
            Lparietal_std = squeeze(std_all(subj, trial, 6, :))';
            Ltemporal_std = squeeze(std_all(subj, trial, 7, :))';
            Loccipital_std = squeeze(std_all(subj, trial, 8, :))';
            
            overall_std = squeeze(mean(std_all(subj, trial, :, :), 3))';
            Rhemi_std = squeeze(mean(std_all(subj, trial, 1:4, :), 3))';
            Lhemi_std = squeeze(mean(std_all(subj, trial, 5:8, :), 3))';
            
            % Add current trial values to spreadsheet
            bands_sheet(currentline,:) = [subj, trial, Rfrontal, Rparietal, Rtemporal, Roccipital, Lfrontal, Lparietal, Ltemporal, Loccipital, Rhemi, Lhemi, overall];
            bands_sheet_std(currentline,:) = [subj, trial, Rfrontal_std, Rparietal_std, Rtemporal_std, Roccipital_std, Lfrontal_std, Lparietal_std, Ltemporal_std, Loccipital_std, Rhemi_std, Lhemi_std, overall_std];
            
            currentline = currentline + 1;
            
        end
        
    end
end

csvwrite('/Volumes/File Storage/JoVE/Data/bands_sheet.csv', bands_sheet);
csvwrite('/Volumes/File Storage/JoVE/Data/bands_sheet_std.csv', bands_sheet_std);

%% Coherence - Isolate only included trials
load('/Volumes/File Storage/JoVE/Data/coherence.mat');
load('/Volumes/File Storage/JoVE/Data/includedlist.mat');

included_ind = 1;
for trials = 1:size(included_list,1)
    cohere_included(included_ind, :,:,:) = squeeze(cohere_allpairs(included_list(trials,1), included_list(trials,2), :,:,:));
    included_ind = included_ind+1;
end

save(['/Volumes/File Storage/JoVE/Data/coherence_included.mat'], 'cohere_included','cohere_freqs', '-v7.3');

%% Coherence - reorganise
% Trials x comparisons x coherence - 111 x 28 x whatever
% (-> use to isolate bands, make 2D spreadsheet, etc.)

% Set current position to beginning of spreadsheet
currentline = 1;

% Cycle through each trial
for subj = 1:size(cohere_included,1)
    
    trial_ind = 1;
    
    for elecone = 1:7
        for electwo = elecone+1:8
            coherence_organised(currentline, trial_ind, :) = squeeze(cohere_included(subj,elecone,electwo,:));
            trial_ind = trial_ind+1;
        end
    end
    
    
    currentline = currentline+1;
    
end


%% Coherence - split into bands


% Theta:
theta_start_sample = find(cohere_freqs==3);
theta_end_sample = find(cohere_freqs==7);

coherence_theta = coherence_organised(:,:,theta_start_sample:theta_end_sample);
avg_coherence_theta = mean(coherence_theta,3);
std_coherence_theta = std(coherence_theta,0,3);

% Alpha:
alpha_start_sample = find(cohere_freqs==7);
alpha_end_sample = find(cohere_freqs==14);

coherence_alpha = coherence_organised(:,:,alpha_start_sample:alpha_end_sample);
avg_coherence_alpha = mean(coherence_alpha,3);
std_coherence_alpha = std(coherence_alpha,0,3);

% Beta:
beta_start_sample = find(cohere_freqs==14);
beta_end_sample = find(cohere_freqs==20);

coherence_beta = coherence_organised(:,:,beta_start_sample:beta_end_sample);
avg_coherence_beta = mean(coherence_beta,3);
std_coherence_beta = std(coherence_beta,0,3);

% Delta:
delta_start_sample = find(cohere_freqs==0.5);
delta_end_sample = find(cohere_freqs==3);

coherence_delta = coherence_organised(:,:,delta_start_sample:delta_end_sample);
avg_coherence_delta = mean(coherence_delta,3);
std_coherence_delta = std(coherence_delta,0,3);

% Gamma:
gamma_start_sample = find(cohere_freqs==20);
gamma_end_sample = find(cohere_freqs==40);  % Changing to 30 as fs = 64 (-> Nyquist)

coherence_gamma = coherence_organised(:,:,gamma_start_sample:gamma_end_sample);
avg_coherence_gamma = mean(coherence_gamma,3);
std_coherence_gamma = std(coherence_gamma,0,3);

% Overall:
avg_coherence_overall = mean(coherence_organised,3);
std_coherence_overall = std(coherence_organised,0,3);

save(['/Volumes/File Storage/JoVE/Data/coherence_organised.mat'], 'coherence_organised', 'coherence_theta', 'coherence_alpha', 'coherence_beta', 'coherence_delta', 'coherence_gamma', 'cohere_freqs', '-v7.3');

%% Coherence - Set up spreadsheets for export

coherence_sheet = [included_list, avg_coherence_theta, avg_coherence_alpha, avg_coherence_beta, avg_coherence_delta, avg_coherence_gamma, avg_coherence_overall];
coherence_std_sheet = [included_list, std_coherence_theta, std_coherence_alpha, std_coherence_beta, std_coherence_delta, std_coherence_gamma, std_coherence_overall];
            
csvwrite('/Volumes/File Storage/JoVE/Data/coherence_sheet.csv', coherence_sheet);
csvwrite('/Volumes/File Storage/JoVE/Data/coherence_std_sheet.csv', coherence_std_sheet);