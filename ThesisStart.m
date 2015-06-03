%%%%%%%%%%%%
%
%   Thesis Work
%
%%%%%%%%%%

%% Add path to BCILAB
addpath(genpath('/home/amart/BCILAB-master'));

% Bring in AAR toolbox
%eeglabRoot = fileparts(which('eeglab'));
%url = 'https://github.com/germangh/eeglab_plugin_aar/archive/master.zip';
%unzip(url, [eeglabRoot filesep 'plugins']);
%addpath(genpath(eeglabRoot));
%eeglab redraw;

%% Start bcilab and eeglab
bcilab;
eeglab;

%cd /home/amart/
%addpath(genpath('/home/amart/BCILAB-master'));
% Clear workspace
clc;

% ThesisTest;