function start_toolkit()
%START_TOOLKIT Loads all the HVACsim Toolkit folder paths in MATLAB. 
%Run this function before calling 'simSystem.m' and the HVACsim modules.
%
% Syntax:  start_toolkit
%
% Inputs:
%    none
%
% Outputs:
%    none
%
% Example: 
%    start_toolkit
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author        : Panayiotis Papadopoulos
% Work address  : KIOS Research Center, University of Cyprus
% email         : papadopoulos.m.panayiotis@gmail.com
% Website       : http://www.kios.ucy.ac.cy
% Last revision : September 2016

%------------- BEGIN CODE --------------
addpath(genpath(pwd));
disp('HVACsim Toolkit Paths Loaded.');    
%------------- END OF CODE --------------

