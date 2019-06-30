# NetworkAnalysis

Scripts for network analysis of EEG data
Allows derivation of network models & analysis based on models

Includes:
NetworkAnalysis_Demonstration.m - walkthrough of data loading through to derivation of models, followed by example analysis including data-driven identification of subgroups, splitting into groups, characterisation of networks within groups, comparison of measures & focused analysis of specific frequency bands

NA_results2spreadsheet.m - converts derived measures into spreadsheets for storing & extracting data easily - adds subject IDs & trial numbers to columns 1 & 2 to allow specific subjects / trials to be easily extracted & compared

NA_reshapecoherence.m - converts coherence measures into matrices for plotting using R scripts

produce_plots.r - R script for producing visualisations of coherence measures & network models using ggplot2 library; also includes several other functions for visualising other parameters
