% This file contains various options used by different MATLAB scripts
% (Flatcone, IMPS, nplot, and others).
% 
% Please copy this file in your active working directory before editing


%   Content:
%   --------
%   - General
%   - Normalization
%   - Data transfer
%   - Averaging, binning and cell size (2D plots)
%   - Plotting
%   - Spectrometer setup
%   - IMPS specific options
%	- Physical constants



%---------------------------
% General:
%---------------------------

% Correction for efficiency of detector channels:
vanacorr = 2;   % 0: none, 
                % 1: take sum of channels in file
                % 2: efficiencies explicitly given below
vanafile = ...  % Filename with vanadium scan
    '062712';
det_eff = ...   % give 31 numbers describing the efficiency of the detector channels
 [1.015,1.1048,0.9938,0.828,1.0554,1.0807,0.95,0.9925,0.9839,0.9408,1.0116,0.9279,0.9593,0.8332,1.0111,0.9753,0.9335,1.0275,0.9984,1.0494,1.0006,1.057,0.9973,1.0213,0.9953,1.026,1.0255,1.0135,1.0124,1.0968,1.082];



% Convention for direction of Q 
% +1 if Q is defined as Q=kf-ki,   (FLATCONE, MAD3d !!)
% -1 for Q=ki-kf                   (MAD !!) 
QSign = +1;


% Use of a3 or a3p...
stdscanmode = 'A3=0';


allchannels = 1:31;
% Which channels of Multidetector to consider?
channels = allchannels;


%---------------------------
% Normalization:
%---------------------------

% Variable on which to normalize (usually 'M1' for Monitor; 'TIME' for time)
normalizeto = 'M1';

% Desired normalization value for this variable (monitor counts, seconds, ...)
normval = 4000;


%---------------------------
% Data transfer
%---------------------------

knownservers = {'in3','in8','in12','thales','in20'};
defaultdirectory = '/users/data';
defaultserver = 'thales';
defaultuser = 'thales';


%---------------------------
% Averaging, grid, binning:
%---------------------------


% Maximum admissible deviations in different variables when combining
% various datasets in which this value should be equal:
maxdeviate.KF = 0.01;
maxdeviate.KI = 0.05;
maxdeviate.QVERT = 0.02;
maxdeviate.TEMP = 1000;
maxdeviate.normal = [.05,.05];
maxdeviate.c = .05;


% Standard grid spacing along the different dimensions for different types
% of coordinates: 
% (During averaging, points are assigned ("binned") to a regular grid to
% overcome insignificant deviations in their coordinate values. This
% defines the "resolution".)
stdgrid.ANGLES = [0.1, 0.1];             % a4'/psi
stdgrid.A4ENERGY = [0.1, .25];          % a4'/en
stdgrid.ANGLESENERGY = [0.1, .25, .1];   %a4'/en/psi
stdgrid.ANGLESQZ = [stdgrid.ANGLES, .02]; %a4'/psi/Qvert
stdgrid.QXY = [0.01, 0.01];             % Qx/Qy
stdgrid.QXQYEN = [stdgrid.QXY, stdgrid.ANGLESENERGY(2)]; % Qx/Qy/En
stdgrid.QXYZ = [stdgrid.QXY, .02];      % Qx/Qy/Qz
stdgrid.LINEARQ = [0.002,0.05];               % Q/en
stdgrid.HKLET = [.005,.005,.01,.01,.005];
stdgrid.SCANSTEP = [.1,.1];
stdgrid.GENERAL = .5;

% Maximum distances for binning
% (often not used because grid fine enough)
% Be careful when using values significantly different from instrumental resolution
stdbindist.ANGLES = [.1, .1];     % a4'/psi
stdbindist.A4ENERGY = [.1, .1];   % a4'/en
stdbindist.ANGLESENERGY = stdbindist.A4ENERGY([1,2,1]);
stdbindist.ANGLESQZ = [stdbindist.ANGLES, .05];
stdbindist.QXY = [0.02, 0.02];    % Qx/Qy
stdbindist.QXQYEN = [stdbindist.QXY, stdbindist.ANGLESENERGY(2)]; % Qx/Qy/En
stdbindist.QXYZ = [stdbindist.QXY, .05];
stdbindist.HKLET = stdgrid.HKLET;
stdbindist.LINEARQ = stdgrid.LINEARQ *5;
stdbindist.SCANSTEP = [.3,.3];
stdbindist.GENERAL = stdgrid.GENERAL;

% For computation of Voronoi cells, how to scale axes in case of different
% dimensions
stdratio.ANGLES = [1,1];     
stdratio.A4ENERGY = [1/30, 1]; 
stdratio.ANGLESENERGY = stdratio.A4ENERGY([1,2,1]);
stdratio.ANGLESQZ = [stdratio.ANGLES, 1];
stdratio.QXY = [1,1]; 
stdratio.QXQYEN = [1,1,1];
stdratio.QXYZ = [1,1,1];
stdratio.LINEARQ = [1, 1/5];
stdratio.QEPLANE = [1, .2];
stdratio.SCANSTEP = [1,1];
stdratio.GENERAL = [1,1];

% Maximum cell sizes along the coordinate axes when setting up the coloured patch for plotting 
% (Using high values fills "unesthetic" gaps, but this may be physically misleading or produce
%  artifacts at the edges)
stdcell.ANGLES = [2.5,2.5];       % a4'/psi
stdcell.A4ENERGY = [2.5,1];      % a4'/en
stdcell.ANGLESENERGY = [stdcell.A4ENERGY, stdcell.ANGLES(2)];
stdcell.ANGLESQZ = [stdcell.ANGLES, .3];
stdcell.QXY = [.16, .16];  stdcell.QPLANE = stdcell.QXY;
stdcell.QXQYEN = [stdcell.QXY, 1];
stdcell.QXYZ = [stdcell.QXY, .3];
stdcell.LINEARQ = [.13, .7];
stdcell.QEPLANE = [.1,1];
stdcell.SCANSTEP = [1.05,1.05];
stdcell.GENERAL = [1,1];
stdcell.HKLVECTORS = [.5,.5];

% ---------------------------
% Plotting:
% ---------------------------

% Title to be displayed in plots etc.
exptitle = 'No Title';

% Standard options
plotopt.showvectors = 1;    %Display orienting vectors?
plotopt.showgrid = 1;       %Display HKL coordinate grid?
plotopt.showedges = 1;      % Show edges of each slice?
plotopt.linlog = 'LIN';     %'lin' for linear, 'log' for log color scale
plotopt.presentation = 'voronoi';  % choose plot type below
%    'voronoi';     patch of Voronoi cells
%    'linear';      linearly (or not) interpolated data
%    'contourf';    Filled contour plot
plotopt.interpolationtype = 'patchinterp';  % Choose interpolation method
%    'pcolorsmooth';    interpolated shading with pcolor grid (ngrid^2 pixels)
%    'pcolorfacet';     flat shading with pcolor grid (ngrid^2 pixels)
%    'patchinterp';     use patch command to shade based on triangulation
%                   of orig. data points (no grid applies; generally looks less nice)
plotopt.interpolationgrid = 200; % Pixels in each dimension if interpolation on a grid 
%                               (large numbers time consuming)
plotopt.interpolationsystem = 'data';   % In which coordinate system to work for interpolation?
%   'data': work on plotstruct.dataset.coordlist;
%   'plot': work on plotstruct.coordlist;
plotopt.interpolationlimit = 2; % Defines when to leave blank a region in the interpolation when real data points too far
%                   (Triangles exceeding this size are thrown out; in units of pixelsize (stdsize)). Typical ~2.
plotopt.interpolationalgorithm = 'natural';
%                   Method used in "griddata" (linear/cubic/natural/nearest)
plotopt.showcells = false;   % Display edges of Voronoi cells?
plotopt.showcoordpoints = false; % Display coordinates of measured points?
plotopt.preferHKL = true;   % HKL coords rather than Angstroms?
plotopt.axisHKL = true;  % In q-plots, use HKL rather than rec.Ang. if possible?



%---------------------------
% Spectrometer setup:
%---------------------------

% Instrumental zeros 
% DO NOT depend on sample alignment, only on instrument!!
% When the encoder offsets are properly set, these values are 0,0,-180 !!
zgu_0 = 0;   % flat goniometer for gu=0
zgl_0 = 0;   % flat goniometer for gl=0
za3_0 = -180;% gl-axis parallel ki for a3=0


% Zeros (values for which the plane defined by ax,bx is horizontal in Lab space with ax||ki)
% DOES depend on sample alignment, used to compensate the "misalignment"
% ** normally replaced by zeros stored in scan files; only for test purposes !! **
stdzeros.gu = 0;
stdzeros.gl = 0;
stdzeros.a3 = -180;
stdzeros.a3p = -180;
% Sample scattering direction (normally replaced by value in scan file)
stdSS = -1;

% Rotation directions for Goniometers for a3=0 with za3 at its nominal value (given by za3_0) 
GUsign = -1; % (ILL-TAS: -1) +1 if gu>0 tilts cryostat top towards reactor for a3=0; -1 otherwise.        IN20 before 2009: +1
GLsign = 1;  % (ILL-TAS: +1) +1 if gl>0 is clockwise rotation around ki seen from reactor; -1 otherwise.  IN20 before 2009: -1

%----------------------------
% powxbu parameters
%----------------------------
powparam.ti = 2;
powparam.da4 = .15;
powparam.np = 30;


%----------------------------
% IMPS
%----------------------------

impsdetcenter = 26.5;  % Central channel (= on a6 arm) of IMPS multidetector (Should normally be 27)
impsanadist = 4; %cm  Distance of IMPS Analyzers (axes)
impsda = 3.266;  %Angstrom, DA for Ge-Analyzer

impsdetupper = 212;
impsdetlower = 46;      % Upper and lower limit of the detector region which is effectively usable

impsroi.stdwidth = 1.6;    % Default width of ROI when automatically calculated
impsroi.mindist = 1;       % Minimum number of "free" channels between two ROIs when automatically calculated (-1 for no check)
impsroi.usefixed = false;  % Use fixed number of channels instead?
impsroi.fixednum = 2;      % Width of ROIs in channels (if used)

%----------------------------
% Directories
%----------------------------
directories.povoutput = 'C:\Dokumente und Einstellungen\Paul\Eigene Dateien\Steffens-Archiv\plots\povray\Dichten\Versuche';

%---------------------------
% Some physical constants
%---------------------------

hbar    = 1.05457168e-34; %J.s
mass_n  = 1.67492428e-27; %kg
meVJ    = 6.24150947e21; %meV/J




