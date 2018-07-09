function test()

clc; clear all; close all;
load('testPER_Temp.mat');
addpath('./Plugins/');
c = physconst('LightSpeed');
fc = 60e9;
lambda = c / fc;

possible_locations = getElementPosition(arrayHandle);
for id = 1:problem.nUsers
    relevant_positions = (W(id,:)~=0);
    Taper_user = W(id,relevant_positions);
    
    handle_Conf_Array_USER{id} = phased.ConformalArray(...
        'Element',arrayHandle.Element,...
        'ElementPosition', [possible_locations(1,relevant_positions);...
        possible_locations(2,relevant_positions);...
        possible_locations(3,relevant_positions)],...
        'Taper',Taper_user);
    resp{id} = phased.ArrayResponse( ...
        'SensorArray', handle_Conf_Array_USER{id}, ...
        'WeightsInputPort', true);
end

%% Pure math
elem_loc = getElementPosition(handle_Conf_Array_USER{1});
angles = [problem.phiUsers; problem.thetaUsers];
elementPattern = cosd(angles(1, 1)) .^ handle_Conf_Array_USER{1}.Element.CosinePower(1) ...
    * cosd(angles(2, 1)) .^ handle_Conf_Array_USER{1}.Element.CosinePower(2);

T_vec = 1 / lambda * elem_loc.' * getSecondMatrixInEq27(angles(1, 1), angles(2, 1)); % Equation 27
v_vec = exp(-1i * 2 * pi * T_vec); % Equation 26
r_from_pure_math = abs(handle_Conf_Array_USER{1}.Taper * v_vec * elementPattern).^2

%% patternAzimuth
r_from_patternAzimuth = patternAzimuth(handle_Conf_Array_USER{1}, ...
             problem.freq,problem.thetaUsers(1),'Azimuth',problem.phiUsers(1),'Type','power')
         
%% array response
resp_from_ArrayRepsonse = abs(resp{1}(fc, angles(:, 1), ones(8, 1))) .^ 2

end

function M = getSecondMatrixInEq27(azang, elang)
M = [-cosd(elang).*cosd(azang);...
    -cosd(elang).*sind(azang);...
    -sind(elang)];
end