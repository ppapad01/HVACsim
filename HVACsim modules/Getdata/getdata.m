function DATA=getdata(name)
%getdata clear DATA.mat and reads the gbXML file from the Benchmark file and returns the mat-file DATA. 
%Run this function after calling 'start_toolkit'.
%
% Syntax:  getdata('name.xml')
%
% Inputs:
%    name='name.xml'
%
% Outputs:
%    DATA : is structure of the gbXML file
%
% Example: 
%    getdata('House_8x23.xml')
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
clear DATA;

DATA=parseXML(name);



