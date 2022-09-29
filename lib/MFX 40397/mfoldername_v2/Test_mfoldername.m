%--------------------------------------------------------------------------
% Test file for mfoldername.m fucntion
% Creates a path to easily save or load your results no matter where your 
% function is located
% Shows three general examples
%--------------------------------------------------------------------------
% Author: Daniel R. Herber, Graduate Student, University of Illinois at
% Urbana-Champaign
% Date: 2/19/2013
%--------------------------------------------------------------------------
% 06/09/2013 v1.1 Documentation changes and example to reflect the functions 
%                 ability to aid with loading data
% 01/04/2016 v1.4 Name change to mfoldername
% 01/11/2016 v2.0 Complete rewrite based on Jan Simon's suggestions
%--------------------------------------------------------------------------
%% Run your code
x = 0;
y = 1;

%% Get path
% 1. You want the current function or script path in new folder
path1 = mfoldername(mfilename('fullpath'),'Saved_Data');
% 2. You want a specific function or script path in function folder
path2 = mfoldername('Test_mfoldername','');
% Extra. You want to load data that is in the function folder
path3 = mfoldername(mfilename('fullpath'),'');

%% Display paths to see the output
disp(path1)
disp(path2)
disp(path3)

%% Use the path to save the results
save([path1,'Data1.mat'],'x','y') % note where Data1.mat is saved
save([path2,'Data2.mat'],'x','y') % note where Data2.mat is saved
disp('Data saved')

%% Use the path to load the data
load([path3,'mfoldername_test_data.mat'])
disp('Data loaded')