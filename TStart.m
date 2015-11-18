%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TStart
%       File for starting EEGLAB
%       Includes changing/adding proper directories (if needed)
%
%		* Currently using EEGLAB version 13.4.4b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
 
%% Start EEGLAB
% If EEGLAB is not already in path it needs to be added or you need to cd into proper directory
eeglab;

 
%% Bring in AAR toolbox - once brought in and redrawn no need to repeat
% eeglabRoot = fileparts(which('eeglab'));
% url = 'https://github.com/germangh/eeglab_plugin_aar/archive/master.zip';
% unzip(url, [eeglabRoot filesep 'plugins']);
% addpath(genpath(eeglabRoot));
% eeglab redraw;

%% Run Proper Script for processing and classification
% TMain 
% TMainWave
TMainWave3
