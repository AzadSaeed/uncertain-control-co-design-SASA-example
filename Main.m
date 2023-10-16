% This function produces results for the following formulations:
% 1. Deterministic CCD, 
% 2. Stochastic in Expectation UCCD,
% 3. Worst-case Robust UCCD 
% for the simple SASA UCCD problem introduced in the following article:
% 
% S. Azad and D. R. Herber, "Investigations into Uncertain Control Co-Design
% Implementations for Stochastic in Expectation and Worst-case Robust "
% in ASME IMECE 2022, Columbus, OH

% For probabilistic uncertaintites, MCS and gPC are used.
% For deterministic representation of uncertainties, uncertainties belong
% to a deterministic set.
% Specific case studies are listed below:
%
% pcase = 0; % Deterministic problem  
%
% pcase = 1; % Stochastic, Nested, Multi-Control, MCS 
% pcase = 2; % Stochastic, Nested, Multi-Control, gPC 
%
% pcase = 3; % SC Worst-case robust
% pcase = 4; % SC Worst-case robust using penalty instead of constraints
%
% pcase = 5; % MC Worst-case- Polytopic uncertainties
% pcase = 6; % MC Worst-case- Uncertainty size (std)
% pcase = 7; % MC Worst-case- Uncertainty size (nsigma)
%
% pcase = 8; % Closed-loop investigations
% pcase = 9; % MPC 
%
% The Simple SASA problem can  be found in the following reference:
% D. R. Herber and J. T. Allison, "Unified Scaling of Dynamic Optimization
% Design Formulations," in Volume 2A: 43rd Design Automation Conference,
% 2017, doi: 10.1115/detc2017-67676
%
% The results from this implementation are published in IMECE 2022.
%
% Contibutor: Saeed Azad, Ph.D.

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


% Runs all cases, saves results and plots figures
clear; clc; close all


% Add the required content
Commonfilesetup;


%-------------------------------------------------------------------------%
% Solve the deterministic problem
% pcase = 0;
% SimpleSASA_UCCD(pcase);
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Solve the Stochastic, Nested, Multi-Control, using MCS  
% pcase = 1;
% SimpleSASA_UCCD(pcase);
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Solve the Stochastic, Nested, Multi-Control, using gPC  
% pcase = 2;
% SimpleSASA_UCCD(pcase);
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Solve the SC, Worst-case robust with partially relaxed Constraints  
% pcase = 3;
% SimpleSASA_UCCD(pcase);
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Solve the SC, Worst-case robust with penalty terms   
% Pen = 0:0.1:1;
% pcase = 4;
% for i = 1:length(Pen)
%       SimpleSASA_UCCD(pcase,Pen(i));
% end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Solve the MC, Worst-case robust using polytopic uncertainties 
% pcase = 5;
% SimpleSASA_UCCD(pcase);
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Solve the MC, Worst-case robust for various sizes of uncertainties
% pcase = 6;
% auxdata.GF = 0:0.05:2;
% for i = 1:length(auxdata.GF)
%     SimpleSASA_UCCD(pcase,[],auxdata.GF(i));
% end

%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Solve the MC, Worst-case robust for various sizes of uncertainties
% pcase = 7;
% auxdata.idx = 0.1:0.1:5;     % constraint shift index
% for i = 1:length(auxdata.idx)
%     SimpleSASA_UCCD(pcase,[],[],auxdata.idx(i));
% end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Closed-loop investigations
% pcase = 8;
% Closed_loop_investigations
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% MPC investigations
% pcase = 9;
% SimpleSASA_UCCD(pcase);

%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Create plots for different cases 0 -7 
% pcase = 7;
% Createplots(pcase)
%-------------------------------------------------------------------------%