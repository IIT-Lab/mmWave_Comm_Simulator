classdef LTE5GPackageInfo < matlab.addons.internal.SupportPackageInfoBase

%   Copyright 2016-2017 The MathWorks, Inc.

    methods 
        function obj = LTE5GPackageInfo()
            %LTE5GPackageInfo Construct package info object
            obj.baseProduct = 'LTE System Toolbox';
            obj.displayName = '5G Library for LTE System Toolbox';
            obj.name = '5G Library for LTE System Toolbox';
            
            sproot = matlabshared.supportpkg.getSupportPackageRoot();
            
            % Files that should always be included in deployed package.
            % Include all support package source files and resources.
            obj.mandatoryIncludeList = {...
                fullfile(sproot, 'toolbox/lte/supportpackages/lte5g/+lte5g') ...
                fullfile(sproot, 'toolbox/lte/supportpackages/lte5g/+nr5g') ...
                fullfile(sproot, 'resources/lte5g')...
                };
            
            % Add conditional includes so that the support package is only
            % suggested if the 5G Library related functions are used.
            obj.conditionalIncludeMap = containers.Map;
            obj.conditionalIncludeMap(fullfile(sproot,'toolbox/lte/supportpackages/lte5g','nr5gTDLChannel.m')) = {};
            obj.conditionalIncludeMap(fullfile(sproot,'toolbox/lte/supportpackages/lte5g','nr5gCDLChannel.m')) = {};
            
        end
        
    end

end