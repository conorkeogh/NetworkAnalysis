%{
EEG Network Analysis

Script for deriving network measures from EEG data

Given list of subjects, loads EEG data, derives network measures & saves
results for analysis

Conor Keogh
Revised 26/06/19
%}

%% Feature Extraction
% Calculate power & coherence measures for each trial

% Define variables for tracking success of analyses
included_line = 1;
excluded_line = 1;

excluded_list=[]; % Prevents crashing if none excluded

for i = 1:length(subjects)

    % Display subject ID
    disp(['Subject ', num2str(i)]);
    
    subj = subjects{i};
    trialnum = trials(i);
    
    for j = 1:trialnum
        
        % Display trial number
        disp(['Trial ', num2str(j)]);
        
        % Load preprocessed data
        data = load([dir 'EEG/Preprocessed/' subj '.mat']);
        eeg = data.eeg_standardised;

        epoch_start = 5 * 60 * 128; % Discarding first five minutes
        epoch_end = 15 * 60 * 128;  % Here using ten minute epoch length
        
        epoch = eeg(:, epoch_start:epoch_end);
        
        try

            epoch = eeg(:, epoch_start:epoch_end);

            % Calculate power spectra spectra:
            [spectra, freqs] = spectopo(epoch, 0, 128);

            % Save spectra:
            all_spectra(i,j,:,:) = spectra;
            
            % Get average power for each frequency band for each channel:

            theta_start_sample = find(freqs==3);
            theta_end_sample = find(freqs==7);

            alpha_start_sample = find(freqs==7);
            alpha_end_sample = find(freqs==14);

            beta_start_sample = find(freqs==14);
            beta_end_sample = find(freqs==20);

            delta_start_sample = find(freqs==0.5);
            delta_end_sample = find(freqs==3);

            gamma_start_sample = find(freqs==20);
            gamma_end_sample = find(freqs==40);

            % Get power in each band for each channel:
            for k = 1:size(spectra, 1)

                disp(['Channel ', num2str(k)]);

                theta_spectrum = spectra(k, theta_start_sample:theta_end_sample);
                [avg_theta, std_theta] = normfit(theta_spectrum);

                alpha_spectrum = spectra(k, alpha_start_sample:alpha_end_sample);
                [avg_alpha, std_alpha] = normfit(alpha_spectrum);

                beta_spectrum = spectra(k, beta_start_sample:beta_end_sample);
                [avg_beta, std_beta] = normfit(beta_spectrum);

                delta_spectrum = spectra(k, delta_start_sample:delta_end_sample);
                [avg_delta, std_delta] = normfit(delta_spectrum);

                gamma_spectrum = spectra(k, gamma_start_sample:gamma_end_sample);
                [avg_gamma, std_gamma] = normfit(gamma_spectrum);
            end
                
                % Get overall:
                [avg_overall, std_overall] = normfit(spectra(k,:));

                % Create data structure containing power measures &
                % standard devidations for each subject, for each trial,
                % for each channel
                powerbands_allsubj(i, j, k, :) = [avg_theta, avg_alpha, avg_beta, avg_delta, avg_gamma, avg_overall];
                std_powerbands_allsubj(i, j, k, :) = [std_theta, std_alpha, std_beta, std_delta, std_gamma, std_overall];
            end
            

        
        % For each unique electrode pair, derive a measure of coherence
        % between the recorded signals
        % Build up a model of statistical relationships between signals
        % across all channels
        for elecone = 1:8
            for electwo = 1:8
                
                [cohere_allpairs(i, j, elecone, electwo, :), cohere_freqs] = mscohere(epoch(elecone, :), epoch(electwo, :), [], [], [], 128);
                
            end 
        end
        
        % List successfully analysed trials:
        included_list(included_line,:) = [i, j];
        included_line = included_line + 1;
            
        catch
            % List trials failed & pad data structure to avoid errors
            disp(['FAILED - ', subj, ' TRIAL ', num2str(j)]);
            powerbands_allsubj(i, j, k, :) = [0, 0, 0, 0, 0];
            std_powerbands_allsubj(i, j, k, :) = [0, 0, 0, 0, 0];
            
            all_spectra(i,j,:,:) = zeros(8,513);
            
            excluded_list(excluded_line,:) = [i, j];
            excluded_line = excluded_line + 1;
        end
    end
end

% Confirm completion in console
disp('DONE');

% Save derived features
avg = powerbands_allsubj;
std_all = std_powerbands_allsubj;
save([dir '/Data/powerbands.mat'], 'avg', 'std_all');
save([dir '/Data/spectra.mat'], 'all_spectra', 'freqs');
save([dir '/Data/coherence.mat'], 'cohere_allpairs', 'cohere_freqs', '-v7.3');
save([dir '/Data/includedlist.mat'], 'included_list', 'excluded_list');
disp('Saved');