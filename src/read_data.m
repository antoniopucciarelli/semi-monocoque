% this program reads the monocoque's data from a txt file
% the data MUST be expressed in this way
%
% POINTS || points #
% || point # || x-coord || y-coord || stringer area ||
% 
% PANELS || panels #
% || panel # || constitutive points || thickness ||
% 
% CELLS  || cells #
% || cell #  || constitutive panels ||
%
% FORCE
% || x-direction || y-direction || Mt || x-coord || y-coord ||
%
% SELECTION
%   selected_point => point which remains out the computation of shear flux
%                         --> because it's not needed
%   
%   selected_panels => panels where the translational equilibrium equations
%                         are not used
%   
% SIGN
% matrices where all the elements represents the contribution of the shear
% flux on the total moment of the structure --> includes q_i and q_star_i
%   GLOBAL-STRUCTURE
%   || PANEL(1) ........................... PANEL(end) || 
%   || +1 => if the induced moment is counterclockwise ||
%   || -1 => if the induced moment is clockwise        ||
%   ith CELL
%   || CELL(ith).PANEL(1) ....... CELL(ith).PANEL(end) || 
%   || +1 => if the induced moment is counterclockwise ||
%   || -1 => if the induced moment is clockwise        ||
%

%% data

% POINTS: 6
POINTSdata = [ ...
    1, 0, 0, 100; ...
    2, 1, 0, 100; ...
    3, 1, 1, 100; ...
    4, 0, 1, 100; ...
    5, 2, 0, 100; ...
    6, 2, 1, 100  ...
];

% PANELS: 6
PANELSdata = [  ...
    1, 1, 2, 0; ...
    2, 2, 3, 0; ...
    3, 3, 4, 0; ...
    4, 4, 1, 0; ... % selected PANEL4
    5, 5, 2, 0; ...
    6, 6, 5, 0; ... % selected PANEL6
    7, 6, 3, 0  ...
];

% CELLS: 2
CELLSdata1 = [1,2,3,4];
CELLSdata2 = [2,5,6,3];

% FORCE
Tx = 0;
Ty = 1800;
Mt = 0;
rx = 1;
ry = 1;

% SELECTION
%   points
    selected_point  = 4;
%   panels
    selected_panels = [4,6];

% % SIGN
% sign_values = [-1,0,-1,0,1,0,-1; ... % GLOBAL 
%                -1,-1,-1,0,0,0,0; ... % LOCAL CELL1
%                 0,0,0,0,1,0,-1,1;... % LOCAL CELL2
% ];
% 
% sign_values_load = [ 0,0,0,-1,0,1,0;    ... % GLOBAL
%                      -1,-1,-1,-1,0,0,0; ... % LOCAL
% ];




