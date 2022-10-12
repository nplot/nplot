function materials = powderdefinitions

% Definition of powder standard samples (used by powxbu).
% Add new materials accordingly.

% P. Steffens, 7/2014


% Silicon
Si.name = 'Silicon';
Si.lattice = [5.4309, 5.4309, 5.4309, 90, 90, 90];
Si.reflections = [1 1 1; 2 2 0; 3 1 1; 4 0 0; 3 3 1; 4 2 2; 5 1 1; 4 4 0; 5 3 1; 6 2 0; 5 3 3; 4 4 4; 5 5 1; 7 1 1; 6 4 2];

% YiG (Y3Fe5O12)
YiG.name = 'YiG';
YiG.lattice = [12.376, 12.376, 12.376, 90, 90, 90];
YiG.reflections = [2 1 1; 2 2 0; 3 2 1; 4 0 0; 4 2 0; 3 3 2; 4 2 2; 4 3 1; 5 2 1; 4 4 0; 6 1 1; 4 4 4; 6 4 0; 5 5 2; 6 4 2; 8 0 0];

% Al2O3
Al2O3.name = 'Al2O3';
Al2O3.lattice = [4.7589, 4.7589, 12.9889, 90, 90, 120];
Al2O3.reflections = [0 1 2; 1 0 4; 0 0 6; 1 1 3; 0 2 4; 1 1 6; 2 1 1; 0 1 8; 2 1 4; 3 0 0; 3 0 6; 2 2 3; 0 0 12; 2 2 6];


%% 
materials = {Si,YiG,Al2O3};

