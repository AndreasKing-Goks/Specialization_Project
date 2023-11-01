function [M_T] = Transform(Matrix, r_o_bg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform()                                                             %
%                                                                         %              
% Transform the a Matrix from it's current frame to the origin frame,     %
% given the distance vector from the current frame to the origin frame.   %
%                                                                         %
% Argument:                                                               %
% Matrix    : Matrix that wanted to be transformed. Dim(6x6).             %  
% r_o_bg    : Distance vector. Dim(3x1).                                  %      
%                                                                         %
% Created:      27.09.2023	Andreas Sitorus                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transformation Matrix
H = H_(r_o_bg);

M_T = H.' * Matrix * H;
end