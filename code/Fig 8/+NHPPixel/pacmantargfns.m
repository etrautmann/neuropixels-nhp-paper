%% PACMANTARGFNS pacman task target functions
% Simulates target functions used in the PacMan task
%
% SYNTAX
%   y = pacmantargfns(Conds, nTrials, varargin)
%
% REQUIRED INPUTS
%   Cond (table) - condition table. Required fields:
%       .type (char): 'STA' (static), 'RMP' (ramp), 'SIN' (sinsoid), 'CHP' (chirp)
%       .offset (scalar) in [0,1]
%       .amplitude (scalar) in [0,1]
%       .duration (scalar) seconds
%       .frequency (2-vector) Hz
%       .nCycles (scalar)
%
%   nTrials (scalar) - number of trials per condition
%
% OPTIONAL INPUTS
%   optIn (class) - description
%
% VARIABLE INPUTS
%   (...,'dt',<scalar>) - time step [sec] (default: 1e-3)
%   (...,'padDur',<scalar>) - pad duration before and after the target
%       (default: 0)
%
% OUTPUTS
%   y (class) - description
%
% EXAMPLE(S) 
%
%   % table headings
%   headings = {'type','offset','amplitude','duration','frequency','nCycles'};
%
%   % condition parameters
%   condParams = {
%       'STA',0.5,NaN,1,NaN,NaN;
%       'RMP',0,1,2,NaN,NaN;
%       'SIN',0,0.5,NaN,1,4;
%       'CHP',0,0.75,6,[0.5,1],NaN
%       };
%
%   % make table
%   C = cell2table(condParams, 'VariableNames', headings);
%
%   % generate targets with 5 trials per condition
%   [y,t] = pacmantargfns(C,5);
%
%
% IMPLEMENTATION
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% SEE ALSO:

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated:

function [y,t,trBounds,condNo, timeFirstPeak] = pacmantargfns(Cond, nTrials, varargin)
% timeFirstPeak is time in seconds where the first peak occurs, as in t time vec.
    
    %% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'PACMANTARGFNS';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'Cond', @istable);
addRequired(P, 'nTrials', @isscalar);
% addOptional(P, 'optIn', default, validationFunction);
addParameter(P, 'dt', 1e-3, @isscalar);
addParameter(P, 'padDur', 0.5, @(x) isscalar(x) || length(x)==2);
addParameter(P, 'iti', 1, @isscalar);
addParameter(P, 'OutputFormat', 'vector', @ischar)
addParameter(P, 'randomize', true, @islogical)

% clear workspace (parser object retains the data while staying small)
parse(P, Cond, nTrials, varargin{:});
clear ans varargin

timeFirstPeak = 1;

%%

% extract parameters
dt = P.Results.dt;
padDur = P.Results.padDur;
if isscalar(padDur)
    padDur = repmat(padDur,1,2);
end

% ensure proper formatting
if ~iscell(Cond.offset)
    Cond.offset = num2cell(Cond.offset);
end
if ~iscell(Cond.amplitude)
    Cond.amplitude = num2cell(Cond.amplitude);
end
if ~iscell(Cond.frequency)
    Cond.frequency = mat2cell(Cond.frequency,ones(height(Cond),1),size(Cond.frequency,2));
end

% loop through condition table
nCond = height(Cond);
y = cell(1,nCond);
for iC = 1:nCond
    
    % determine duration of target function
    targDur = round(Cond.duration(iC),3);
    
    % construct time vector
    t = -padDur(1):dt:(targDur+padDur(2));
    
    % indices of target and post-pad zones
    targIdx = t>=0 & t<=targDur;
    postPadIdx = t>targDur;
    
    % construct target function
    y{iC} = zeros(size(t));
    switch Cond.type{iC}
        case 'RMP'
            A = Cond.amplitude{iC}(1);
            m = A/targDur;
            
            y{iC} = targIdx .* (m*t) + ...
                postPadIdx .* (A+y{iC});

            timeFirstPeak = 1; % arbitrary
            
        case 'TRI'
            A = Cond.amplitude{iC}(1);
            m = A/(targDur/2);
            
            y{iC} = (t>=0 & t<=(targDur/2)) .* (m*t) + ...
                (t>(targDur/2) & t<=targDur) .* (-m*t + m*targDur);
            
        case 'POW'
            A = Cond.amplitude{iC}(1);
            pow = Cond.power;
            
            y{iC} = A*(t/t(end)).^pow;
            
        case 'SIN'
            
            nSin = length(Cond.frequency{iC});
            
            for jS = 1:nSin
                om = 2*pi*Cond.frequency{iC}(jS);
                A = Cond.amplitude{iC}(jS);
                
                y{iC} = y{iC} + ...
                    targIdx .* (A/2 * (1-cos(om*t))) + ...
                    postPadIdx .* (A/2 * (1-cos(om*targDur))); % + y{iC});
            end

            % first peak where om*t == pi/2 --> t = pi /2 / om
            om = 2*pi*Cond.frequency{iC}(1);
            assert(nnz(Cond.frequency{iC}) == 1, "multiple freq not implemented yet")
            timeFirstPeak = t(find(targIdx .* om.*t >= pi, 1, 'first'));
            
        case 'CHP'
            k = diff(Cond.frequency{iC}(1:2))/targDur;
            f0 = Cond.frequency{iC}(1);
            A = Cond.amplitude{iC}(1);
            iFin = floor(1+(padDur(1)+targDur)/dt); % changed per Eric's email
            
            y{iC} = targIdx .* (A/2 * (1-cos(2*pi*t.*(f0+k/2*t)))) + ...
                postPadIdx .* ((A/2 * (1-cos(2*pi*t(iFin)*(f0+k/2*t(iFin))))) + y{iC});

            timeFirstPeak = t(find(targIdx .* (2*pi*t.*(f0+k/2*t)) >= pi, 1, 'first'));
    end
    
    y{iC} = y{iC} + sum(Cond.offset{iC});
end

if P.Results.randomize
    % shuffle condition numbers by trial
    nCond = height(Cond);
    condNo = cell2mat(cellfun(@randperm, num2cell(nCond*ones(1,nTrials)), 'uni', false));
    
    % call each target function by its condition number
    y = y(condNo);
end

condNo = condNo';

if strcmp(P.Results.OutputFormat, 'vector')
    % intersperse an ITI between each trial
    y = [y; repmat({zeros(1,P.Results.iti/dt)}, 1, length(y))];
    
    % reshape and remove last ITI
    y = reshape(y,1,numel(y));
    y(end) = [];
    
    % trial bounds
    trBounds = [0 cumsum(cellfun(@length,y))];
    trBounds = reshape(trBounds, 2, length(trBounds)/2)';
    
    y = cell2mat(y);
    
    % generate time vector
    t = (1:length(y))*dt - dt;
    
    % orient vertically
    y = y';
    t = t';
    
elseif strcmp(P.Results.OutputFormat, 'cell')
    y = y(:);
    t = cellfun(@(x) (1:length(x))*dt - (dt+padDur(1)), y, 'uni', false);
    
    % orient vertically
    y = cellfun(@(x) x(:), y, 'uni', false);
    t = cellfun(@(x) x(:), t, 'uni', false);
    
    trBounds = [];
    
else
    y = [];
    t = [];
end