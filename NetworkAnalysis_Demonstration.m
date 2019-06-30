%{
EEG Network Analysis

Script for demonstrating loading EEG data, deriving network model,
visualising networks & applying analysis techniques to models

Some discrete functions, available in provided library, are expanded here
for demonstration purposes

Following demonstration:
Loads & prepares data
Derives measures
Uses measures to identify groups within data
Splits data into groups
Analyses & compares networks between identified groups
Performs focused analysis on overall spectrum

Conor Keogh
Revised 26/06/19
%}

%% Prepare to load data

% Set directory containing EEG data
dir = '/Volumes/File Storage/JoVE/';

% Load list of EEG files to be analysed
file = [dir '/Data/file_list.csv'];
fid = fopen(file);
trial_details = textscan(fid, '%s %f', 'Delimiter', ',');  % Edit to read in any other details in file as required, e.g. groups, trial number, etc.
subjects = trial_details{1}; % Derive any further variables from trial_details structure as required for grouping / analysis
trials = ones(length(subjects),1); % If variable numbers of trials for each subject, derive from file list - here just using one trial per subject

%% Data Cleaning
%{
Visually inspect data, perform artefact rejection & epoch data. Re-examine
data to ensure appropriate quality. This should be a semi-manual procedure,
as described in the accompanying text.

For the purposes of demonstration, assume that data has been cleaned &
epoched appropriately.
%}

%% Preprocess each EEG file

% Load each EEG data in file_list & preprocess
for i = 1:length(subjects)
    
    % Isolate current subject
    subj = subjects{i};
    disp(['Processing subject: ' subj]);
    
    % Load EEG data
    file = [dir '/EEG/' subj '.csv'];
    eeg = csvread(file);
    eeg = eeg'; % Data may need to be transposed, depending on the format it is saved in
    
    % Perform preprocessing; function available as NA_preprocess.m
    % Implementation as function: eeg = NA_preprocess(eeg); 
    data_shifted=[];
    data_rerefed=[];
    data_filtered=[];
    
    % Correct baseline for each channel:
    for i = 1:size(eeg,1)
        baseline = mean(eeg(i,:),2);
        data_shifted(i,:) = eeg(i,:) - baseline;
    end
    
    % Rereference to average:
    data_avg = mean(data_shifted,1);
    
    for i = 1:size(eeg,1)
        data_rerefed(i,:) = data_shifted(i,:) - data_avg;
    end
    
    % Apply filters:
    % High pass:
    [H_HP, G_HP] = butter(4, (2*1)/128, 'high');
    % Low pass:
    [H_LP, G_LP] = butter(4, (2*40)/128, 'low');
    
    for i = 1:size(eeg,1)
        data_filtered(i,:) = filtfilt(H_HP, G_HP, data_rerefed(i,:));
        data_filtered(i,:) = filtfilt(H_LP, G_LP, data_filtered(i,:));
    end
    
    eeg = data_filtered;
    
    % Save preprocessed data
    save([dir '/EEG/Preprocessed/' subj '.mat'], 'eeg');
    disp('Saved');
    
end

% Confirm completion in console
disp('DONE');

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

        % Define start & end; will depend on earlier epoching
        epoch_start = 1;
        epoch_end = 10 * 60 * 128;  % Here using ten minutes
        try

            epoch = eeg(:, epoch_start:epoch_end);

            % Calculate power spectra spectra: (requires EEGLab loaded)
            [spectra, freqs] = spectopo(epoch, 0, 128);

            % Save spectra:
            all_spectra(i,j,:,:) = spectra;
            
            % Get average power for each frequency band for each channel:

            % Theta:
            theta_start_sample = find(freqs==3);
            theta_end_sample = find(freqs==7);

            % Alpha:
            alpha_start_sample = find(freqs==7);
            alpha_end_sample = find(freqs==14);

            % Beta:
            beta_start_sample = find(freqs==14);
            beta_end_sample = find(freqs==20);

            % Delta:
            delta_start_sample = find(freqs==0.5);
            delta_end_sample = find(freqs==3);

            % Gamma:
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
                % Replace w pwelch & cpsd
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

%% Convert data structures into spreadsheets to facilitate analysis

%{
In order to facilitate analysis, we recommend converting derived features
into spreadsheets; this allows easy visualisation of data to ensure no
errors have occurred during analysis & allows data to be easily extracted
for grouping by condition, trial, etc. for comparisons

A script demonstrating how to restructure these data into easily usable
spreadsheets is included in NA_results2spreadsheet.m; the details are
excluded here for demonstration purposes

Note that it is possible to perform comparisons between groups etc. without
the intermediate step of converting to spreadsheets, but analysis is
rendered more straightforward by performing the conversion
%}

%% Clustering
% Derive distance measure within feature space
% Use measure to identify groups within data
% Use these groups to inform further analyses

eucD = pdist(all_coherence(:,3:end), 'euclidean');   % Derive distance metric
clustTree = linkage(eucD, 'ward');

[h,nodes] = dendrogram(clustTree,0);
h_gca = gca;
h_gca.TickDir = 'out';
% Appears to cluster into 2 main groups

% Split subjects into two groups based on dendrogram:
cluster_inds = cluster(clustTree, 'maxclust', 2); % Cluster indices for each subject
group1 = find(cluster_inds == 1);
group2 = find(cluster_inds == 2);

%% Compare groups
% Load extracted features, divide into groups & compare based on models

% Load data
all_power = csvread([dir '/Data/bands_sheet.csv']);
all_coherence = csvread([dir '/Data/coherence_sheet.csv']);

% Set up groups: extract model data based on subject group if known groups
% Commented out here to allow grouping based on clustering
%group1 = [1,2,6,7,9,10]; % Example groupings; edit based on whatever way data being analysed is grouped
%group2 = [3,4,5,8];

g1_power = all_power(group1, 3:end); % Extract all data related to this group; exclude column 1 & 2, which correspond to subject ID and trial ID in spreadsheet
g1_coherence = all_coherence(group1, 3:end);

g2_power = all_power(group2, 3:end);
g2_coherence = all_coherence(group2, 3:end);

% Any individual measures of interest can be directly compared at this
% stage using standard statistical tests
% Performing individual comparisons of all model parameters results in a
% large number of comparisons, so is not well suited to anything but
% exploratory analyses
% An example of a systematic exploratory analysis is given in
% NA_comparegroups.m

% Visualise coherence patterns:
% Isolate coherence measures within each band & reshape for visualisation
g1_coherence_matrix = mean(g1_coherence,1);
g1_coherence_matrix_theta = g1_coherence_matrix(1:28);
g1_coherence_matrix_alpha = g1_coherence_matrix(29:56);
g1_coherence_matrix_beta = g1_coherence_matrix(57:84);
g1_coherence_matrix_delta = g1_coherence_matrix(65:112);
g1_coherence_matrix_gamma = g1_coherence_matrix(113:140);
g1_coherence_matrix_overall = g1_coherence_matrix(141:168);

g2_coherence_matrix = mean(g2_coherence,1);
g2_coherence_matrix_theta = g2_coherence_matrix(1:28);
g2_coherence_matrix_alpha = g2_coherence_matrix(29:56);
g2_coherence_matrix_beta = g2_coherence_matrix(57:84);
g2_coherence_matrix_delta = g2_coherence_matrix(65:112);
g2_coherence_matrix_gamma = g2_coherence_matrix(113:140);
g2_coherence_matrix_overall = g2_coherence_matrix(141:168);

% Function NA_reshapecoherence maps coherence data onto a grid for
% visualisation
g1_coherence_matrix_theta = NA_reshapecoherence(g1_coherence_matrix_theta);
g1_coherence_matrix_alpha = NA_reshapecoherence(g1_coherence_matrix_alpha);
g1_coherence_matrix_beta = NA_reshapecoherence(g1_coherence_matrix_beta);
g1_coherence_matrix_delta = NA_reshapecoherence(g1_coherence_matrix_delta);
g1_coherence_matrix_gamma = NA_reshapecoherence(g1_coherence_matrix_gamma);
g1_coherence_matrix_overall = NA_reshapecoherence(g1_coherence_matrix_overall);

g2_coherence_matrix_theta = NA_reshapecoherence(g2_coherence_matrix_theta);
g2_coherence_matrix_alpha = NA_reshapecoherence(g2_coherence_matrix_alpha);
g2_coherence_matrix_beta = NA_reshapecoherence(g2_coherence_matrix_beta);
g2_coherence_matrix_delta = NA_reshapecoherence(g2_coherence_matrix_delta);
g2_coherence_matrix_gamma = NA_reshapecoherence(g2_coherence_matrix_gamma);
g2_coherence_matrix_overall = NA_reshapecoherence(g2_coherence_matrix_overall);

% Save data for plotting
csvwrite([dir '/Data/g1_coherence_matrix_theta.csv'], g1_coherence_matrix_theta);
csvwrite([dir '/Data/g1_coherence_matrix_alpha.csv'], g1_coherence_matrix_alpha);
csvwrite([dir '/Data/g1_coherence_matrix_beta.csv'], g1_coherence_matrix_beta);
csvwrite([dir '/Data/g1_coherence_matrix_delta.csv'], g1_coherence_matrix_delta);
csvwrite([dir '/Data/g1_coherence_matrix_gamma.csv'], g1_coherence_matrix_gamma);
csvwrite([dir '/Data/g1_coherence_matrix_overall.csv'], g1_coherence_matrix_overall);

csvwrite([dir '/Data/g2_coherence_matrix_theta.csv'], g2_coherence_matrix_theta);
csvwrite([dir '/Data/g2_coherence_matrix_alpha.csv'], g2_coherence_matrix_alpha);
csvwrite([dir '/Data/g2_coherence_matrix_beta.csv'], g2_coherence_matrix_beta);
csvwrite([dir '/Data/g2_coherence_matrix_delta.csv'], g2_coherence_matrix_delta);
csvwrite([dir '/Data/g2_coherence_matrix_gamma.csv'], g2_coherence_matrix_gamma);
csvwrite([dir '/Data/g2_coherence_matrix_overall.csv'], g2_coherence_matrix_overall);

%% Data Visualisation
%{
Coherence data can be visualised directly in MATLAB; however, the included
R scripts produce better visualisations

It is further possible to perform statistical analyses on all model
parameters & visualise all those with a statistically significant result
for exploratory analysis. This is demonstrated in NA_comparegroups.m
%}

%% Network Comparisons
% Calculate covariance measures & visualise networks
% Use covariance measures to derive eigenvectors for PCA
% Perform DR & compare PCs

% Get covariances matrices for all & groups:
g1_cov = cov(g1_coherence); % Group 1 coherence - for visualisation
g2_cov = cov(g2_coherence); % Group 2 coherence - for visualisation

all_pca = pca(all_coherence(:,3:end)');  % Get PCA coefficients for all subjects
% Note that PCA typically requires a large sample size for stable results;
% for the present demonstration we are using a small sample to limit
% computation time
g1_pca = all_pca(group1,:); % Isolate group 1
g2_pca = all_pca(group2,:); % Isolate group 2

% Plot groups on new axes based on principal components:
figure; hold on;
plot(g1_pca(:,1), g1_pca(:,2), 'b.', 'MarkerSize', 30);
plot(g2_pca(:,1), g2_pca(:,2), 'r.', 'MarkerSize', 30);
legend('Group 1', 'Group 2');
xlabel('PC 1');
ylabel('PC 2');

% Compare networks based on first principal component:
p = ranksum(g1_pca(:,1), g2_pca(:,1)); % Compare first principal components
disp(['Group 1 vs. Group 2: P = ' num2str(p) ', Wilcoxon rank sum test']);

% Network covariances can then be visualised using included R scripts
% Or alternatively:
network_columns = horzcat('Theta', repmat({' '},1,27), 'Alpha', repmat({' '},1,27), 'Beta', repmat({' '},1,27), 'Delta', repmat({' '},1,27), 'Gamma', repmat({' '},1,27), 'Overall', repmat({' '},1,27));
HeatMap(g1_cov, 'Colormap', 'redbluecmap', 'RowLabels', network_columns, 'ColumnLabels', network_columns, 'ColumnLabelsLocation', 'top', 'ColumnLabelsRotate', 45);
HeatMap(g2_cov, 'Colormap', 'redbluecmap', 'RowLabels', network_columns, 'ColumnLabels', network_columns, 'ColumnLabelsLocation', 'top', 'ColumnLabelsRotate', 45);

csvwrite([dir '/Data/g1_cov.csv'], g1_cov);
csvwrite([dir '/Data/g2_cov.csv'], g2_cov);
%% ROI Analysis
% Isolate data for region of interest
% Derive network visualisations & perform DR for isolated data

% Select just data within overall spectrum:
all_coherence_overall = all_coherence(:,143:end);   % Can select any anatomic / frequency criteria for selecting data
g1_coherence_overall = all_coherence_overall(group1,:);
g2_coherence_overall = all_coherence_overall(group2,:);

% Perform analyses as above using isolated data:
% Get covariances matrices for all & groups:
g1_cov_overall = cov(g1_coherence_overall); % Group 1 coherence - for visualisation
g2_cov_overall = cov(g2_coherence_overall); % Group 2 coherence - for visualisation

all_pca_overall = pca(all_coherence_overall');  % Get PCA coefficients for all subjects

g1_pca_overall = all_pca_overall(group1,:); % Isolate group 1
g2_pca_overall = all_pca_overall(group2,:); % Isolate group 2

% Plot groups on new axes based on principal components:
figure; hold on;
plot(g1_pca_overall(:,1), g1_pca_overall(:,2), 'b.', 'MarkerSize', 30);
plot(g2_pca_overall(:,1), g2_pca_overall(:,2), 'r.', 'MarkerSize', 30);
legend('Group 1', 'Group 2');
xlabel('PC 1');
ylabel('PC 2');

% Compare networks based on first principal component:
p = ranksum(g1_pca_overall(:,1), g2_pca_overall(:,1)); % Compare first principal components
disp(['Group 1 vs. Group 2, Overall Spectrum Only: P = ' num2str(p) ', Wilcoxon rank sum test']);

% Network covariances can then be visualised using included R scripts
% Or alternatively:
covariancelabels = {'Fp2-C4', 'Fp2-O2', 'Fp2-T4', 'Fp2-Fp1', 'Fp2-C3', 'Fp2-T3', 'Fp2-O1', ...
                          'C4-O2', 'C4-T4', 'C4-Fp1', 'C4-C3', 'C4-T3', 'C4-O1', 'O2-T4', ...
                          'O2-Fp1', 'O2-C3', 'O2-T3', 'O2-O1', 'T4-Fp1', 'T4-C3', 'T4-T3', ...
                          'T4-O1', 'Fp1-C3', 'Fp1-T3', 'Fp1-O1', 'C3-T3', 'C3-O1', 'T3-O1'};
HeatMap(g1_cov_overall, 'Colormap', 'redbluecmap','RowLabels', covariancelabels, 'ColumnLabels', covariancelabels,'ColumnLabelsLocation', 'top', 'ColumnLabelsRotate', 90);
HeatMap(g2_cov_overall, 'Colormap', 'redbluecmap','RowLabels', covariancelabels, 'ColumnLabels', covariancelabels,'ColumnLabelsLocation', 'top', 'ColumnLabelsRotate', 90);

csvwrite([dir '/Data/g1_cov_overall.csv'], g1_cov_overall);
csvwrite([dir '/Data/g2_cov_overall.csv'], g2_cov_overall);
