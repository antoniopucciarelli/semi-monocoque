% this program computes the stress distribution on beam-like structure using the 
% semimonocoque method 
%
% hp: 
% 	-- panel thickness ~ 0 --> panel shear stress constant over panel thickness
%   -- normal traction applied only on the stringers --> panel shear stress variation null over the panel length
%   -- the panels shear stress can change only passing the stringers 
%

clc
clear 
close all

% generating structure
%  |-> points generation
%  |-> panels generation
%  |-> cells generation
%  |-> 1st order inertia moment 
%  --> 2nd order inertia moment 

%% points 

% opening file where are stored all the data needed
FILE = fopen('semimonocoque.txt','r');

% scanning file for POINT's data
if (fscanf(FILE, '%s', 1) == "POINTS")
    
    Npoints = fscanf(FILE, '%d', 1);      % # of POINTS inside the structure
    POINTS  = POINTSobj.empty(Npoints,0); % POINT object allocation
    
    for ii = 1:Npoints
        POINTS(ii).coords(1) = fscanf(FILE, '%f', 1);
        POINTS(ii).coords(2) = fscanf(FILE, '%f', 1);
        POINTS(ii).area      = fscanf(FILE, '%f', 1);
    end
    
else
    error("invalid txt file --> POINTS");
end

%% panels 

% number of panels
if(fscanf(FILE, '%s', 1) == "PANELS")
    
    Npanels = fscanf(FILE, '%d', 1);      % # of panels inside the structure
    PANELS  = PANELSobj.empty(Npanels,0); % PANEL object allocation
    
    % allocating panel constitutive points and thickness data
    for ii = 1:Npanels
        PANELS(ii).points = fscanf(FILE, '%f', 2);
        PANELS(ii).thick  = fscanf(FILE, '%f', 1);
        PANELS(ii)        = PANELS(ii).get_panels(POINTS); % panels' points allocation
        PANELS(ii)        = PANELS(ii).compute_len();      % panels' lenght
    end
    
end

% allocating POINTS object's PANELS belongings
for ii = 1:Npoints
    
    % PANELS' START point  --> the shear flux enters inside the pole
    % PANELS' FINISH point --> the shear flux exits from the pole
    % [ counter_in, counter_out ] = vector dimension counter --> in order to avoid warnings
    % on the increasing size of the position_in and position_out vector
    
    position_in  = zeros(Npanels, 1); % allocation of space for the points number
    position_out = zeros(Npanels, 1); % allocation of space for the points number
    
    counter_in   = 1;
    counter_out  = 1;
    
    for jj = 1:Npanels
        if(PANELS(jj).points(1) == ii)
            position_in(counter_in) = jj;
            counter_in              = counter_in + 1;
        end
        if(PANELS(jj).points(2) == ii)
            position_out(counter_out) = jj;
            counter_out               = counter_out + 1;
        end
    end
    
    POINTS(ii).panels_in  = position_in(1:counter_in-1);
    POINTS(ii).panels_out = position_out(1:counter_out-1);

end

% POINTS' area computation
for ii = 1:Npoints
    for jj = 1:length(POINTS(ii).panels_in)
        POINTS(ii).area = POINTS(ii).area + PANELS(POINTS(ii).panels_in(jj)).len / 2 * PANELS(POINTS(ii).panels_in(jj)).thick;
    end
    for jj = 1:length(POINTS(ii).panels_out)
        POINTS(ii).area = POINTS(ii).area + PANELS(POINTS(ii).panels_out(jj)).len / 2 * PANELS(POINTS(ii).panels_out(jj)).thick;
    end
end

%% cells generation 
if (fscanf(FILE, '%s', 1) == "CELLS")
    
    Ncells = fscanf(FILE, '%d', 1);    % # of cells in the structure
    CELLS  = CELLSobj.empty(Ncells,0); % generating cells object array 
    
    for ii = 1:Ncells
        SIZEcell         = fscanf(FILE, '%d', 1);        % # of panels that forms the cell
        CELLS(ii).panels = fscanf(FILE, '%d', SIZEcell); % cell's panel
    end
    
else
    error("invalid txt file --> CELLS");
end

%% force and moment loading 
if (fscanf(FILE, '%s', 1) == "LOAD")
    
    Tx = fscanf(FILE, '%f', 1);
    Ty = fscanf(FILE, '%f', 1);
    Mt = fscanf(FILE, '%f', 1);
    rx = fscanf(FILE, '%f', 1);
    ry = fscanf(FILE, '%f', 1);
    
else
    error("invalid txt file --> LOAD");
end
%% computing structural geometry properties

% computing baricenter
CG   = [0,0];
Atot = 0;

for ii = 1:Npoints
	CG   = CG   + POINTS(ii).CG_comp();
	Atot = Atot + POINTS(ii).area;
end

CG = CG / Atot; 

% translation of points wrt the baricenter
for ii = 1:Npoints
    POINTS(ii) = POINTS(ii).change_coords(CG);
end
	
% computing inertia coefficients
% 1st order 
for ii = 1:Npoints
	POINTS(ii) = POINTS(ii).compute_Sx();
	POINTS(ii) = POINTS(ii).compute_Sy(); 
end

% 2nd order
Ixx = 0;
Iyy = 0;
Ixy = 0;
 
for ii = 1:Npoints
	Ixx = Ixx + POINTS(ii).area * POINTS(ii).coords(2)^2;
	Iyy = Iyy + POINTS(ii).area * POINTS(ii).coords(1)^2;
	Ixy = Ixy + POINTS(ii).area * POINTS(ii).coords(1) * POINTS(ii).coords(2);
end

% initializing rotation matrix
ROT = [1, 0; ...
       0, 1];

% computing 2nd time the Sx, Sy, Ixx, Iyy and Ixy values for the structure	 
if(abs(Ixy) > 1e-15)
    
    % computing alpha angle 
    alpha = 0.5 * atan(2*Ixy/(Iyy - Ixx));
    
    % rotating points coords 
    ROT = [ cos(alpha), sin(alpha); ...
           -sin(alpha), cos(alpha)];

    for ii = 1:length(POINTS)
        POINTS(ii).coords = ROT * POINTS(ii).coords';
    end

    % computing inertia coefficients
    % 1st order 
    for ii = 1:length(POINTS)
        POINTS(ii).Sx = POINTS(ii).area * POINTS(ii).coords(2);
        POINTS(ii).Sy = POINTS(ii).area * POINTS(ii).coords(1);
    end

    % 2nd order
    Ixx = 0;
    Iyy = 0;
    Ixy = 0;

    for ii = 1:length(POINTS) 
        Ixx = Ixx + POINTS(ii).area * POINTS(ii).coords(2)^2;
        Iyy = Iyy + POINTS(ii).area * POINTS(ii).coords(1)^2;
        Ixy = Ixy + POINTS(ii).area * POINTS(ii).coords(1) * POINTS(ii).coords(2);
    end
    
else 
    
    alpha = 0;

end

% POINTS rotation
for ii = 1:Npoints
   POINTS(ii).coords = ROT * POINTS(ii).coords'; 
end

% getting new points' coordinates on the panel
for ii = 1:Npanels
      PANELS(ii) = PANELS(ii).get_panels(POINTS);
end

% modifies on load
LOAD = ROT * [Tx; Ty];
Tx   = LOAD(1);
Ty   = LOAD(2);

% modofies on load application
rx = rx - CG(1);
ry = ry - CG(2);
R  = ROT * [rx; ry];
rx = R(1);
ry = R(2);

% printing inertia coefficients
fprintf('Ixx   = %d\nIyy   = %d\nIxy   = %d\nalpha = %d\nCG    = [%f, %f]\n\n', Ixx, Iyy, Ixy, alpha/pi*180, CG(1), CG(2));

% with the change of reference system --> CG == [0; 0]
CG = zeros(2);

%% PLOT - data
plt_structure(POINTS,PANELS,CG,Npoints,Npanels);

%% loading structural analysis data
if (fscanf(FILE, '%s', 1) == "SELECTED_PANELS")
    
    Nsel_panels     = fscanf(FILE, '%d', 1);
    selected_panels = fscanf(FILE, '%d', Nsel_panels);
else 
    error("invalid txt file --> SELECTED PANELS");
end

if (fscanf(FILE, '%s', 1) == "SELECTED_POINT")    
    selected_point = fscanf(FILE, '%d', 1);
else   
    error("invalid txt file --> SELECTED POINT");
end

%% solution for the translation problem
% sum(q) = - Ty * Sx/Ixx * y - Tx * Sy/Iyy * x

% plotting FORCE
figure(1)
S = 1e+4; % scaling force vector in rappresentation
quiver(rx,ry,Tx/S,Ty/S,'-b','LineWidth',3);

MATRIX = zeros(length(POINTS),length(PANELS));
FORCE  = zeros(length(POINTS),1);

% matrix assembling procedure
for ii = 1:length(POINTS)
    for jj = 1:length(POINTS(ii).panels_in)
        MATRIX(ii,POINTS(ii).panels_in(jj))  = -1;
    end
    for jj = 1:length(POINTS(ii).panels_out)
        MATRIX(ii,POINTS(ii).panels_out(jj)) = 1;
    end
end

% known vector assembling procedure 
for ii = 1:length(POINTS)
    FORCE(ii) = - Ty * POINTS(ii).Sx / Ixx - Tx * POINTS(ii).Sy / Iyy;
end

% elimination of the redundant row/s for the computation of the translation
% equilibrium
MATRIX(:,selected_panels) = [];
MATRIX(selected_point,:)  = [];
FORCE(selected_point)     = [];

q = MATRIX \ FORCE;

% allocating fluxes in panel array objects
% fluxes --> panel_in  >> 1st point of PANELS' points --> takes -1 in MATRIX
%            panel_out >> 2nd point of PANELS' points --> takes  1 in MATRIX
% shear stress flux starts from the end point of the panel and ends 
%   at the starting point of the panel
for ii = 1:length(PANELS)
    PANELS(ii).flux = 0;
end

counter = 0;
for ii = 1:length(q)
    if(ii ~= selected_panels)
        PANELS(ii+counter).flux = q(ii);
    else
        PANELS(ii+counter).flux = 0; 
        counter = counter + 1;
    end
end

fprintf('translation problem --> #q = %d\n', length(q));
for ii = 1:length(PANELS)
    fprintf('q%d = %d\n', ii, PANELS(ii).flux);
end

%% computing circulating shear fluxes
M = Ty * rx + Tx * (- ry) + Mt;

fprintf('\nmoment generated from shear forces and torque\nM = %d\n\nLOADS\nTx       = %d\nTy       = %d\nMt       = %d\nT-coords = [%d, %d]\n', M, Tx, Ty, Mt, rx, ry);

if Ncells > 1
    
    % collecting data
    Qmatr = zeros(Ncells * (Ncells+1), Npanels); % matrix for the allocation of the coefficients

    % q* coefficients --> [CELLS(1:Ncells); GLOBAL] study
    % 1:Ncells ==> CELLS  data
    % Ncells+1 ==> GLOBAL data
    for ii = 1:Ncells

        qname = append('q', num2str(ii), '*');

        if(fscanf(FILE, '%s', 1) == qname)
            for jj = 1:Ncells+1
                Qmatr((ii-1)*(Ncells+1) + jj,:) = fscanf(FILE, '%d', Npanels);
            end
        else 
            error("invalid txt file --> " + qname);
        end

    end

    % collecting data for the cells shear flux coefficients
    for ii = 1:Ncells

        cell_name = append('CELL', num2str(ii));

        if(fscanf(FILE, '%s', 1) == cell_name)

            CELLS(ii).coeffs = fscanf(FILE, '%d', Npanels);

        else 
            error("invalid txt file --> " + cell_name);
        end

    end

    % collecting data for the GLOBAL equilibrium study
    if(fscanf(FILE, '%s', 1) == "Mtot")
        Mvec = fscanf(FILE, '%d', Npanels);
        fclose(FILE);                       % all data collected --> closing file
    else
        error("invalid txt file --> Mtot");    
    end

    % generation of the rotational and moment equilibrium matrix
    MATRIX = zeros(Ncells + 1);
    Vknown = zeros(Ncells + 1, 1);

    % cells rotation study
    for ii = 1:Ncells

        for jj = 1:Ncells
            for kk = 1:Npanels
                MATRIX(ii,jj) = MATRIX(ii,jj) + Qmatr((ii-1)*(Ncells+1)+jj, kk) * PANELS(kk).len;
            end
        end

        for kk = 1:Npanels
           Vknown(ii) = Vknown(ii) - PANELS(kk).flux * PANELS(kk).len * CELLS(ii).coeffs(kk); 
        end

    end

    % global equilibrium
    for jj = 1:Ncells

        for kk = 1:Npanels

            X = [CG(1), PANELS(kk).start.coords(1), PANELS(kk).finish.coords(1)];
            Y = [CG(2), PANELS(kk).start.coords(2), PANELS(kk).finish.coords(2)];

            MATRIX(Ncells+1,jj) = MATRIX(Ncells+1,jj) + 2 * Qmatr(jj*(Ncells+1),kk) * polyarea(X,Y);

        end

    end

    for kk = 1:Npanels

        X = [CG(1), PANELS(kk).start.coords(1), PANELS(kk).finish.coords(1)];
        Y = [CG(2), PANELS(kk).start.coords(2), PANELS(kk).finish.coords(2)];

        Vknown(Ncells+1) = Vknown(Ncells+1) + 2 * Mvec(kk) * PANELS(kk).flux * polyarea(X,Y);

    end

    MATRIX(1:Ncells,Ncells+1) = 1;
    Vknown(Ncells+1)          = - Vknown(Ncells+1) + M;

    q_star = MATRIX \ Vknown;
    
    fprintf('\ncirculating shear stress flux --> #%d\n', Ncells);
    for ii = 1:Ncells
        fprintf('q*(%d) = %d\n', ii, q_star(ii));
    end

    fprintf('\ntheta_dot = %d\n', q_star(Ncells+1));

else 
    
    % collecting data
    Qmatr = zeros(1, Npanels); % matrix for the allocation of the coefficients

    % q* coefficients --> [GLOBAL] study
    % GLOBAL data
    if(fscanf(FILE, '%s', 1) == "q1*")
        Qmatr(1,:) = fscanf(FILE, '%d', Npanels);
    else 
        error("invalid txt file --> " + qname);
    end
    
    % getting cell data
    if(fscanf(FILE, '%s', 1) == "CELL1")
        CELLS(1).coeffs = fscanf(FILE, '%d', length(CELLS(1).panels));
    else 
        error("invalid txt file --> CELL1");
    end
    
    % collecting data for the GLOBAL equilibrium study
    if(fscanf(FILE, '%s', 1) == "Mtot")
        Mvec = fscanf(FILE, '%d', Npanels);
        fclose(FILE);                       % all data collected --> closing file
    else
        error("invalid txt file --> Mtot");    
    end
    
    % initializing data
    MATRIX = 0;
    Vknown = 0;
    
    % global equilibrium MATRIX
    for kk = 1:Npanels
        X = [CG(1), PANELS(kk).start.coords(1), PANELS(kk).finish.coords(1)];
        Y = [CG(2), PANELS(kk).start.coords(2), PANELS(kk).finish.coords(2)];
        MATRIX = MATRIX + 2 * Qmatr(1,kk) * polyarea(X,Y);
    end
    % globel known vector 
    for kk = 1:Npanels
        X = [CG(1), PANELS(kk).start.coords(1), PANELS(kk).finish.coords(1)];
        Y = [CG(2), PANELS(kk).start.coords(2), PANELS(kk).finish.coords(2)];
        Vknown = Vknown + 2 * Mvec(kk) * PANELS(kk).flux * polyarea(X,Y);
    end
    
    Vknown  = - Vknown + M;
    
    q_star = Vknown / MATRIX; 
    
    fprintf('\ncirculating shear stress flux --> #%d\n', Ncells);
    fprintf('q*(%d) = %d\n', 1, q_star);
    fprintf('no theta_dot --> single cell structure\n');
    
    fprintf('\nq in structure\n');
    for ii = 1:Npanels
        PANELS(ii).flux_p = q_star; 
        fprintf('q_tot(%d) = %d\n', ii, PANELS(ii).flux + PANELS(ii).flux_p);
    end
    
end
