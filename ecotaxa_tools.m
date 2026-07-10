% ecotaxa_tools.m
% Entry point: loads EcoTaxa TSV exports, selects the imaging instrument,
% builds the "base" database, and launches spectral processing.

close all
clear all
warning('off', 'all');

[directoryoftoolbox, ~, ~] = fileparts(mfilename('fullpath'));

%% Metadata reminder
f = warndlg('Reminder: ecotaxa tsv export must be "by sample", with all metadata correctly filled in and properly formatted — otherwise the pipeline may fail or give wrong results.', 'Before you start');
uiwait(f);

%% Select folder containing TSV files
f = msgbox('Please select the folder containing the *.tsv files (one file per sample)');
uiwait(f);
folder = uigetdir;
cd(folder);

A = dir('*.tsv');
if isempty(A)
    msgbox('The folder contains no *.tsv files. Please export your data from EcoTaxa with the option "separate by sample".');
    return
end

%% Select instrument
list = {'Zooscan', 'Flowcam', 'IFCB', 'PlanktoScope'};
[indx2, tf] = listdlg('PromptString', 'Which instrument was used?', ...
    'ListString', list, 'SelectionMode', 'single', 'ListSize', [300, 150]);
if ~tf, return; end

%% Load TSV files and build spectral database
switch indx2
    case 1  % ZooScan
        [base, m, zoo_groups, Idlistshort] = process_tsv_zooscan;
        process_spectra_and_save;

    case 2  % FlowCam
        [base, m, zoo_groups, Idlistshort] = process_tsv_flowcam;
        process_spectra_and_save;

    case 3  % IFCB
        [base, m, zoo_groups, Idlistshort] = process_tsv_ifcb;
        process_spectra_and_save_IFCB;

    case 4  % PlanktoScope
        [base, m, zoo_groups, Idlistshort] = process_tsv_planktoscope;
        process_spectra_and_save;
end
