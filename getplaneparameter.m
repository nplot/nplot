function [planevectors, origin] = getplaneparameter(planenormals, planec)

% Return parameter form of subspace defined by intersection of n hyperplanes

% P. Steffens, 03/2009

planevectors = null(planenormals);

origin = linsolve(planenormals,planec(:));