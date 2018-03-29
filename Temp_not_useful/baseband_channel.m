classdef baseband_channel < matlab.System
    % untitled2 Add summary here
    %
    % NOTE: When renaming the class name untitled2, the file name
    % and constructor name must be updated to use the class name.
    %
    % This template includes most, but not all, possible properties,
    % attributes, and methods that you can implement for a System object.
    
    % Public, tunable properties
    properties
        
    end
    
    % Public, non-tunable properties
    properties(Nontunable)
        %% AWGN configurations
        noise_method = 'Signal to noise ratio (SNR)';
        SNR = 14;
        
    end
    
    properties(DiscreteState)
        
    end
    
    % Pre-computed constants
    properties(Access = private)
        awgnChannel;
    end
    
    methods
        % Constructor
        function obj = baseband_channel(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:})
        end
    end
    
    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.awgnChannel = comm.AWGNChannel( ...
                'NoiseMethod', obj.noise_method, ...
                'SNR', obj.SNR);
            
            gc.CyclicPrefixLength = 16;
            gc.Doppler = 5;
            gc.chanSRate = 1e6;
            DelaySpread = gc.CyclicPrefixLength -1;
            numPaths= 1;
            gc.PathDelays = floor(linspace(0,DelaySpread,numPaths))*(1/gc.chanSRate);
            gc.PathGains  = zeros(size(gc.PathDelays));
            for n=2:numPaths
                gc.PathGains(n) = gc.PathGains(n-1)-abs(randn);
            end
            
            obj.MIMOChannel = comm.MIMOChannel(...
                'SampleRate',                   obj.sample_rate, ...
                'MaximumDopplerShift',          gc.Doppler, ...
                'PathDelays',                   floor(linspace(0,DelaySpread,numPaths))*(1/gc.chanSRate), ...
                'AveragePathGains',             gc.PathGains,...
                'NumTransmitAntennas',          obj.Nt, ...
                'NumReceiveAntennas',           obj.Nr, ...
                'PathGainsOutputPort',          true,...
                'NormalizePathGains',           true,...
                'NormalizeChannelOutputs',      true, ...
                'SpatialCorrelationSpecification', 'None');
            %             'TransmitCorrelationMatrix',    eye(obj.Nt),...
            %                 'ReceiveCorrelationMatrix',     eye(obj.Nr),...
            
            lightSpeed = physconst('light');
            obj.lambda = lightSpeed / obj.center_frequency;
        end
        
        function rxBits = stepImpl(obj, txBits)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            release(obj.awgnChannel);
            rxBits = obj.awgnChannel(txBits);
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            release(obj.awgnChannel);
        end
        
        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj
            
            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);
            
            % Set private and protected properties
            %s.myproperty = obj.myproperty;
        end
        
        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s
            
            % Set private and protected properties
            % obj.myproperty = s.myproperty;
            
            % Set public properties and states
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end
        
        %% Advanced functions
        function validateInputsImpl(obj,u)
            % Validate inputs to the step method at initialization
        end
        
        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values
        end
        
        function ds = getDiscreteStateImpl(obj)
            % Return structure of properties with DiscreteState attribute
            ds = struct([]);
        end
        
        function processTunedPropertiesImpl(obj)
            % Perform actions when tunable properties change
            % between calls to the System object
        end
        
        function flag = isInputSizeMutableImpl(obj,index)
            % Return false if input size cannot change
            % between calls to the System object
            flag = false;
        end
        
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object
            % configuration, for the command line and System block dialog
            flag = false;
        end
    end
end
