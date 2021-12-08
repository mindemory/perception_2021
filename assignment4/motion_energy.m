%%  motion energy
% Given oddFast, oddSlow, evenFast, evenSlow, the function computes Even1,
% Odd1, Even2, Odd2, where 1 refers to left/up and 2 refers to right/down
% depending on the orientation of the grating and the direction in which it
% is moving. These Even's and Odd's are then used to compute the
% corresponding Energies by taking a sum of element-wise square of the
% elements.
%%

function [Even1, Odd1, Even2, Odd2, Energy1, Energy2]...
    = motion_energy(oddFast, oddSlow, evenFast, evenSlow)  
    Even1 = oddFast + evenSlow;
    Odd1 = -oddSlow + evenFast;
    Even2 = -oddFast + evenSlow;
    Odd2 = oddSlow + evenFast;
    Energy1 = Even1.^2 + Odd1.^2; 
    Energy2 = Even2.^2. + Odd2.^2;
end