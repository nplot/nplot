function gfc=gfcneu(a4alt,gfcalt,a4neu)

% Find new gfc value which corresponds to the setting (a4alt,gfcalt), i.e.
% same Qvert, when a4 is set to a4neu


gamma = pi/180*gfcalt(:);
delta = pi/180*(a4alt(:)-52.5);
deltaneu = pi/180*(a4neu(:)-52.5);

gammaneu = asin ( cos(delta).*sin(gamma)./cos(deltaneu));

gfc = 180/pi*gammaneu;
