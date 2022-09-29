%--------------------------------------------------------------------------
% DTQP_M.m
% Create sequences for Mayer terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [I,J,V] = DTQP_M(Mfull,in,opts)

% extract some of the variables
nt = in.nt; ini = in.i; I_stored = in.I_stored;

% initialize storage arrays
Isav = {}; Jsav = {}; Vsav = {};

% go through each Mayer term
for k = 1:length(Mfull)

    % obtain current substructure
    Mleft = Mfull(k).left;
    Mright = Mfull(k).right;
    Mmatrix = Mfull(k).matrix;

	% obtain matrix
    Mt = DTQP_tmultiprod(Mmatrix,[],0);

    % check if both left and right fields are equal to 0
    if (Mleft ~= 0), R = ini{Mleft}; else, R = 0; end
    if (Mright ~= 0), C = ini{Mright}; else, C = 0; end

    % determine locations and matrix values at this points
    for i = 1:length(R) % number of row continuous variables
        for j = 1:length(C) % number of column continuous variables

            % get current matrix value
            Mv = Mt(:,i,j);

            % check if this entry is always 0
            if any(Mv)

                % Hessian index sequence
                r = DTQP_getQPIndex(R(i),Mleft,0,nt,I_stored); % row
                c = DTQP_getQPIndex(C(j),Mright,0,nt,I_stored); % column

                % assign main diagonal
                Isav{end+1} = r; % row
                Jsav{end+1} = c; % column
                Vsav{end+1} = Mv; % value

            end

        end % end C for loop
    end % end R for loop
end % end M for loop

% combine
I = vertcat(Isav{:});
J = vertcat(Jsav{:});
V = vertcat(Vsav{:});

end % end DTQP_M function