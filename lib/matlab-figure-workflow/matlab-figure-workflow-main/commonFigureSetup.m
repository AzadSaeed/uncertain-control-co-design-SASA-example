% initialize
silentflag = 0; % silent

% add contents to path
RunSilent('AddSubmissionContents(mfilename)',silentflag)

% download required web zips
RunSilent('RequiredWebZips',silentflag)

% add contents to path (files have been downloaded)
RunSilent('AddSubmissionContents(mfilename)',silentflag)

%--------------------------------------------------------------------------
function RequiredWebFiles

disp('-> Obtaining required web files')

% initialize index
ind = 0;

% initialize structure
files = struct('url','','folder','');

% file 1
ind = ind + 1; % increment
files(ind).url = 'http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/legslbdiff.m';
files(ind).folder = 'legslbdiff';

% file 2
ind = ind + 1; % increment
files(ind).url = 'http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/lepoly.m';
files(ind).folder = 'lepoly';

% obtain full function path
full_fun_path = which(mfilename('fullpath'));
outputdir = fullfile(fileparts(full_fun_path),'lib');

% download
DownloadWebFiles(files,outputdir)

disp(' ')
end
%--------------------------------------------------------------------------
function RequiredWebZips

disp('-> Obtaining required web zips')

% initialize index
ind = 0;

% initialize structure
zips = struct('url','','folder','','test','');

% zip
ind = ind + 1; % increment
zips(ind).url = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/40397/versions/7/download/zip/mfoldername_v2.zip';
zips(ind).folder = 'MFX 40397';
zips(ind).test = 'mfoldername';

% zip
ind = ind + 1; % increment
zips(ind).url = 'https://github.com/altmany/export_fig/archive/master.zip';
zips(ind).folder = 'MFX 23629';
zips(ind).test = 'export_fig';

% obtain full function path
full_fun_path = which(mfilename('fullpath'));
outputdir = fullfile(fileparts(full_fun_path),'lib');

% download and unzip
DownloadWebZips(zips,outputdir)

disp(' ')
end

%--------------------------------------------------------------------------
function AddSubmissionContents(name)

disp('-> Adding submission contents to path')
disp(' ')

% turn off potential warning
warning('off','MATLAB:dispatcher:nameConflict')

% current file
fullfuncdir = which(name);

% current folder
submissiondir = fullfile(fileparts(fullfuncdir));

% add folders and subfolders to path
addpath(genpath(submissiondir))

% turn warning back on
warning('on','MATLAB:dispatcher:nameConflict')
end

%--------------------------------------------------------------------------
function DownloadWebZips(zips,outputdir)

% store the current directory
olddir = pwd;

% create a folder for outputdir
if ~exist(outputdir, 'dir')
    mkdir(outputdir); % create the folder
else
    addpath(genpath(outputdir)); % add folders and subfolders to path
end

% change to the output directory
cd(outputdir)

% go through each zip
for k = 1:length(zips)

    % get data
    url = zips(k).url;
    folder = zips(k).folder;
    test = zips(k).test;

    % first check if the test file is in the path
    if exist(test,'file') == 0

        try
            % download zip file
            zipname = websave(folder,url);

            % save location
            outputdirname = fullfile(outputdir,folder);

            % create a folder utilizing name as the foldername name
            if ~exist(outputdirname, 'dir')
                mkdir(outputdirname);
            end

            % unzip the zip
            unzip(zipname,outputdirname);

            % delete the zip file
            delete([folder,'.zip'])

            % output to the command window
            disp(['Downloaded and unzipped ',folder])

        catch % failed to download
            % output to the command window
            disp(['Failed to download ',folder])

            % remove the html file
            delete([folder,'.html'])
        end

    else
        % output to the command window
        disp(['Already available ',folder])
    end
end

% change back to the original directory
cd(olddir)
end
%--------------------------------------------------------------------------
function RunSilent(str,silentflag)

% if silent, capture the output
if silentflag
    O = evalc(str); %#ok<NASGU>
else
    eval(str);
end
end