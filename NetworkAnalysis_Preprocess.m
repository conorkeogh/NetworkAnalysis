%{
EEG Network Analysis

Script for proprocessing data for analysis

Includes:
Baseline correction
Referencing (to mean of all channels)
Filtering

Returns preprocessed data

Conor Keogh
Revised 26/06/19
%}

function processed_data = NetworkAnalysis_Preprocess(data)

%% Preprocess
% Remove baseline:
for i = 1:size(data,1)
    baseline = mean(data(i,:),2);
    data_shifted(i,:) = data(i,:) - baseline;
end

% Rereference to average:
data_avg = mean(data_shifted,1);

for i = 1:size(data,1)
    data_rerefed(i,:) = data(i,:) - data_avg;
end

% Filter:
% High pass:
[H_HP, G_HP] = butter(4, (2*1)/128, 'high');
% Low pass:
[H_LP, G_LP] = butter(4, (2*50)/128, 'low');

for i = 1:size(data,1)
    data_filtered(i,:) = filtfilt(H_HP, G_HP, data_rerefed(i,:));
    data_filtered(i,:) = filtfilt(H_LP, G_LP, data_filtered(i,:));
end

processed_data = data_filtered;