classdef (StrictDefaults)nr5gTDLChannel < matlab.System
%nr5gTDLChannel TR 38.901 Tapped Delay Line (TDL) channel
%   CHAN = nr5gTDLChannel creates a TDL MIMO fading channel System object,
%   CHAN. This object filters a real or complex input signal through the
%   TDL MIMO channel to obtain the channel impaired signal. This object
%   implements the following aspects of TR 38.901:
%   * Section 7.7.2 Tapped Delay Line (TDL) models
%   * Section 7.7.3 Scaling of delays 
%   * Section 7.7.6 K-factor for LOS channel models
%   * Section 7.7.5.2 TDL extension: Applying a correlation matrix
%
%   CHAN = nr5gTDLChannel(Name,Value) creates a TDL MIMO channel object,
%   CHAN, with the specified property Name set to the specified Value. You
%   can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(CHAN,X) filters input signal X through a TDL MIMO fading
%   channel and returns the result in Y. The input X can be a double
%   precision data type scalar, vector, or 2D matrix with real or complex
%   values. X is of size Ns-by-Nt, where Ns is the number of samples and Nt
%   is the number of transmit antennas. Y is the output signal of size
%   Ns-by-Nr, where Nr is the number of receive antennas. Y is of double
%   precision data type with complex values.
% 
%   [Y,PATHGAINS] = step(CHAN,X) returns the MIMO channel path gains of the
%   underlying fading process in PATHGAINS. This syntax applies when you
%   set the PathGainsOutputPort property of CHAN to true. PATHGAINS is of
%   size Ns-by-Np-by-Nt-by-Nr, where Np is the number of paths, i.e., the
%   length of the PathDelays property value of CHAN. PATHGAINS is of double
%   precision data type with complex values.
% 
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(obj,x) and y = obj(x) are
%   equivalent.
%
%   nr5gTDLChannel methods:
%
%   step           - Filter input signal through a TDL MIMO fading channel
%                    (see above)
%   release        - Allow property value and input characteristics changes
%   clone          - Create TDL channel object with same property values
%   isLocked       - Locked status (logical)
%   reset          - Reset states of filters, and random stream if the
%                    RandomStream property is set to 'mt19937ar with seed'
%   info           - Return characteristic information about the TDL 
%                    channel
%   getPathFilters - Get filter impulse reponses for the filters which 
%                    apply the path delays to the input waveform
%
%   nr5gTDLChannel properties:
%
%   DelayProfile               - TDL delay profile
%   PathDelays                 - Discrete path delay vector (s)
%   AveragePathGains           - Average path gain vector (dB)
%   FadingDistribution         - Rayleigh or Rician fading
%   KFactorFirstTap            - K-factor of first tap (dB)
%   DelaySpread                - Desired delay spread (s)
%   MaximumDopplerShift        - Maximum Doppler shift (Hz)
%   KFactorScaling             - Enable K-factor scaling (logical)
%   KFactor                    - Desired Rician K-factor (dB)
%   SampleRate                 - Input signal sample rate (Hz)
%   MIMOCorrelation            - Correlation between UT and BS antennas
%   Polarization               - Antenna polarization arrangement
%   TransmissionDirection      - Transmission direction (Uplink/Downlink)
%   NumTransmitAntennas        - Number of transmit antennas
%   NumReceiveAntennas         - Number of receive antennas
%   TransmitCorrelationMatrix  - Transmit spatial correlation matrix (or 3-D array)
%   ReceiveCorrelationMatrix   - Receive spatial correlation matrix (or 3-D array)
%   TransmitPolarizationAngles - Transmit polarization slant angles in degrees
%   ReceivePolarizationAngles  - Receive polarization slant angles in degrees
%   XPR                        - Cross polarization power ratio (dB)
%   SpatialCorrelationMatrix   - Combined correlation matrix (or 3-D array)
%   NormalizePathGains         - Normalize path gains (logical)
%   PathGainsOutputPort        - Enable path gain output (logical)
%   InitialTime                - Start time of fading process (s)
%   NumSinusoids               - Number of sinusoids in sum-of-sinusoids technique
%   RandomStream               - Source of random number stream
%   Seed                       - Initial seed of mt19937ar random number stream
%   NormalizeChannelOutputs    - Normalize channel outputs (logical)
% 
%   % Example 1: 
%   %   Create an LTE waveform for Reference Measurement Channel (RMC) R.9
%   %   (20MHz, 64QAM, R=0.75), pass it through a channel with delay 
%   %   profile TDL-C, 300ns delay spread and UT velocity 30km/h, and plot
%   %   the received waveform spectrum:
%   
%   rmc = lteRMCDL('R.9');
%   rmc.TotSubframes = 1;
%   data = [1; 0; 0; 1];
%   [txWaveform,~,txInfo] = lteRMCDLTool(rmc,data);
%   
%   v = 30.0;                    % UT velocity in km/h
%   fc = 4e9;                    % carrier frequency in Hz
%   c = physconst('lightspeed'); % speed of light in m/s
%   fd = (v*1000/3600)/c*fc;     % UT max Doppler frequency in Hz
%
%   tdl = nr5gTDLChannel;
%   tdl.NumTransmitAntennas = size(txWaveform,2);
%   tdl.DelayProfile = 'TDL-C';
%   tdl.DelaySpread = 300e-9;
%   tdl.MaximumDopplerShift = fd;
%   tdl.SampleRate = txInfo.SamplingRate;
%   tdl.NumReceiveAntennas = 2;
%
%   rxWaveform = tdl(txWaveform);
%
%   analyzer = dsp.SpectrumAnalyzer('SampleRate',tdl.SampleRate);
%   analyzer.Title = ['Received signal spectrum for ' tdl.DelayProfile];
%   analyzer(rxWaveform);
%
%   % Example 2:
%   %   Plot the path gains of a TDL-E delay profile in a SISO case for a
%   %   Doppler shift of 70Hz:
%    
%   tdl = nr5gTDLChannel;
%   tdl.SampleRate = 500e3;
%   tdl.NumTransmitAntennas = 1;
%   tdl.NumReceiveAntennas = 1;
%   tdl.MaximumDopplerShift = 70;
%   tdl.PathGainsOutputPort = true;
%   tdl.DelayProfile = 'TDL-E';
%    
%   % dummy input signal, its length determines the number of path gain
%   % samples generated
%   in = zeros(1000,tdl.NumTransmitAntennas);
%    
%   % generate path gains
%   [~, pathGains] = tdl(in);
%   mesh(10*log10(abs(pathGains)));
%   view(26,17); xlabel('channel path');
%   ylabel('sample (time)'); zlabel('magnitude (dB)');
%
%   % Example 3:
%   %   Create an LTE RMC R.12 (1.4MHz, QPSK, R=1/3, 4 CRS ports) waveform 
%   %   and pass it through a channel with delay profile TDL-D, 10ns delay 
%   %   spread and a desired overall K-factor of 7.0dB. The channel is
%   %   configured for cross-polar antennas according to TS 36.101 Annnex
%   %   B.2.3A.3 4x2 high correlation.
%
%   rmc = lteRMCDL('R.12');
%   rmc.TotSubframes = 1;
%   data = [1; 0; 0; 1];
%   [txWaveform,~,txInfo] = lteRMCDLTool(rmc,data);
%
%   tdl = nr5gTDLChannel;
%   tdl.NumTransmitAntennas = size(txWaveform,2);
%   tdl.DelayProfile = 'TDL-D';
%   tdl.DelaySpread = 10e-9;
%   tdl.KFactorScaling = true;
%   tdl.KFactor = 7.0; % desired model K-factor (K_desired) dB
%   tdl.SampleRate = txInfo.SamplingRate;
%   tdl.MIMOCorrelation = 'High';
%   tdl.Polarization = 'Cross-Polar';
%
%   rxWaveform = tdl(txWaveform);
%
%   % Example 4:
%   %   Create an LTE RMC R.9 (20MHz, 64QAM, R=0.75) waveform and pass it 
%   %   through a channel with a customized delay profile having two taps 
%   %   as follows: 
%   %   tap 1: Rician, average power 0dB, K-factor 10dB, delay zero
%   %   tap 2: Rayleigh, average power -5dB, delay 45ns
%
%   rmc = lteRMCDL('R.9');
%   rmc.TotSubframes = 1;
%   data = [1; 0; 0; 1];
%   [txWaveform,~,txInfo] = lteRMCDLTool(rmc,data);
%
%   tdl = nr5gTDLChannel;
%   tdl.NumTransmitAntennas = size(txWaveform,2);
%   tdl.DelayProfile = 'Custom';
%   tdl.FadingDistribution = 'Rician';
%   tdl.KFactorFirstTap = 10.0; % K-factor of 1st tap (K_1) in dB
%   tdl.PathDelays = [0.0 45e-9];
%   tdl.AveragePathGains = [0.0 -5.0];
%   tdl.SampleRate = txInfo.SamplingRate;
%
%   rxWaveform = tdl(txWaveform);
%   
%   See also nr5gCDLChannel, lteFadingChannel, comm.MIMOChannel.

%   Copyright 2016-2018 The MathWorks, Inc.
    
%#codegen

% =========================================================================
%   public interface

    methods (Access = public)
        
        % nr5gTDLChannel constructor
        function obj = nr5gTDLChannel(varargin)
            
            % Set property values from any name-value pairs input to the
            % constructor
            setProperties(obj,nargin,varargin{:});
            
        end
        
        function h = getPathFilters(obj)
        %getPathFilters Get path filter impulse reponses
        %   H = getPathFilters(obj) returns a double precision real matrix
        %   of size Nh-by-Np where Nh is the number of impulse response
        %   samples and Np is the number of paths. Each column of H
        %   contains the filter impulse response for each path of the delay
        %   profile. This information facilitates reconstruction of a
        %   perfect channel estimate when used in conjunction with the
        %   PATHGAINS output of the step method (available when the
        %   PathGainsOutputPort property is set to true).
            
            h = nr5gTDLChannel.makePathFilters(obj);            
            
        end
        
    end
    
    % public properties 
    properties (Access = public, Nontunable)
        
        %DelayProfile TDL delay profile
        %   Specify the TDL delay profile as one of 'TDL-A', 'TDL-B',
        %   'TDL-C', 'TDL-D', 'TDL-E' or 'Custom'. See TR 38.901 Section
        %   7.7.2, Tables 7.7.2-1 to 7.7.2-5. When you set this property to
        %   'Custom', the delay profile is configured using the PathDelays,
        %   AveragePathGains, FadingDistribution and KFactorFirstTap
        %   properties.
        %
        %   The default value of this property is 'TDL-A'.
        DelayProfile = 'TDL-A';
        
        %PathDelays Discrete path delay vector (s)
        %   Specify the delays of the discrete paths in seconds as a
        %   double-precision, real, scalar or row vector. This property
        %   applies when DelayProfile is set to 'Custom'.
        %
        %   The default value of this property is 0.0. 
        PathDelays = 0.0;
        
        %AveragePathGains Average path gain vector (dB)
        %   Specify the average gains of the discrete paths in deciBels as
        %   a double-precision, real, scalar or row vector.
        %   AveragePathGains must have the same size as PathDelays. This
        %   property applies when DelayProfile is set to 'Custom'.
        %
        %   The default value of this property is 0.0. 
        AveragePathGains = 0.0;

        %FadingDistribution Fading process statistical distribution
        %   Specify the fading distribution of the channel as one of
        %   'Rayleigh' or 'Rician'. This property applies when DelayProfile
        %   is set to 'Custom'.
        %
        %   The default value of this property is 'Rayleigh', i.e. the
        %   channel is Rayleigh fading. 
        FadingDistribution = 'Rayleigh';
        
        %KFactorFirstTap K-factor of first tap (dB)
        %   Specify the K-factor of the first tap of the delay profile in
        %   dB (K_1) as a scalar. This property applies when DelayProfile
        %   is set to 'Custom' and FadingDistribution is set to 'Rician'.
        %
        %   The default value of this property is 13.3dB. This is the value
        %   defined for delay profile TDL-D.
        KFactorFirstTap = 13.3;
        
        %DelaySpread Desired RMS delay spread in seconds
        %   Specify the desired RMS delay spread in seconds (DS_desired) as
        %   a scalar. See TR 38.901 Section 7.7.3, and Tables 7.7.3-1 and
        %   7.7.3-2 for examples of desired RMS delay spreads. This
        %   property applies when you set the DelayProfile property to
        %   'TDL-A', 'TDL-B', 'TDL-C', 'TDL-D' or 'TDL-E'.
        %
        %   The default value of this property is 30e-9.
        DelaySpread = 30e-9;
        
        %MaximumDopplerShift Maximum Doppler shift (Hz)
        %   Specify the maximum Doppler shift for all channel paths in
        %   Hertz as a double precision, real, nonnegative scalar. The
        %   Doppler shift applies to all the paths of the channel. When you
        %   set the MaximumDopplerShift to 0, the channel remains static
        %   for the entire input. You can use the reset method to generate
        %   a new channel realization.
        %
        %   The default value of this property is 5Hz.
        MaximumDopplerShift = 5.0;
        
    end
        
    properties (Access = public, Nontunable, Logical)
        
        %KFactorScaling Apply K-factor scaling (logical)
        %   Set this property to true to apply K-factor scaling as
        %   described in TR 38.901 Section 7.7.6. Note that K-factor
        %   scaling modifies both the path delays and path powers. This
        %   property applies if DelayProfile is set to 'TDL-D' or 'TDL-E'.
        %
        %   The default value of this property is false.
        KFactorScaling = false;
        
    end
    
    properties (Access = public, Nontunable)
        
        %KFactor Desired Rician K-factor (dB)
        %   Specify the desired K-factor in dB (K_desired) as a scalar.
        %   This property applies when you set the KFactorScaling property
        %   to true. See TR 38.901 Section 7.7.6, and see Table 7.5-6 for
        %   typical K-factors. Note that K-factor scaling modifies both the
        %   path delays and path powers. Note that the K-factor applies to
        %   the overall delay profile i.e. the K-factor after the scaling
        %   is K_model described in TR 38.901 Section 7.7.6, the ratio of
        %   the power of the LOS part of the first path to the total power
        %   of all the Rayleigh paths, including the Rayleigh part of the
        %   first path.
        %
        %   The default value of this property is 9.0dB.
        KFactor = 9.0;
        
        %SampleRate Sample rate (Hz)
        %   Specify the sample rate of the input signal in Hz as a double
        %   precision, real, positive scalar.
        %
        %   The default value of this property is 30.72e6Hz.
        SampleRate = 30.72e6;
        
        %MIMOCorrelation Correlation between UT and BS antennas
        %   Specify the desired MIMO correlation as one of 'Low', 'Medium',
        %   'Medium-A', 'UplinkMedium', 'High' or 'Custom'. Other than
        %   'Custom', the values correspond to MIMO correlation levels
        %   defined in TS 36.101 and TS 36.104. The 'Low' and 'High'
        %   correlation levels are the same for both uplink and downlink
        %   and are therefore applicable to both TS 36.101 and TS 36.104.
        %   Note that 'Low' correlation is equivalent to no correlation
        %   between antennas. The 'Medium' and 'Medium-A' correlation
        %   levels are defined in TS 36.101 Annex B.2.3.2 for
        %   TransmissionDirection = 'Downlink'. The 'Medium' correlation
        %   level is defined in TS 36.104 Annex B.5.2 for
        %   TransmissionDirection = 'Uplink'. When you set this property to
        %   'Custom', the correlation between UE antennas is specified
        %   using the ReceiveCorrelationMatrix property and the correlation
        %   between BS antennas is specified using the
        %   TransmitCorrelationMatrix property. See TR 38.901 Section
        %   7.7.5.2.
        %
        %   The default value of this property is 'Low'.
        MIMOCorrelation = 'Low';
        
        %Polarization Antenna polarization arrangement
        %   Specify the antenna polarization arrangement as one of
        %   'Co-Polar', 'Cross-Polar' or 'Custom'. 
        %
        %   The default value of this property is 'Co-Polar'.
        Polarization = 'Co-Polar';
        
        %TransmissionDirection Transmission direction (Uplink/Downlink)
        %   Specify the transmission direction as one of 'Downlink' |
        %   'Uplink'. This property applies when you set the
        %   MIMOCorrelation property to 'Low', 'Medium', 'Medium-A',
        %   'UplinkMedium' or 'High'.
        %
        %   The default value of this property is 'Downlink'.
        TransmissionDirection = 'Downlink';
        
        %NumTransmitAntennas Number of transmit antennas
        %   Specify the number of transmit antennas as a numeric, real,
        %   positive integer scalar. This property applies when you set the
        %   MIMOCorrelation property to 'Low', 'Medium', 'Medium-A',
        %   'UplinkMedium' or 'High', or when both the MIMOCorrelation and
        %   Polarization properties are set to 'Custom'.
        %
        %   The default value of this property is 1.
        NumTransmitAntennas = 1;
        
        %NumReceiveAntennas Number of receive antennas
        %   Specify the number of receive antennas as a numeric, real,
        %   positive integer scalar. This property applies when you set the
        %   MIMOCorrelation property to 'Low', 'Medium', 'Medium-A',
        %   'UplinkMedium' or 'High'.
        %
        %   The default value of this property is 2.
        NumReceiveAntennas = 2;
        
        %TransmitCorrelationMatrix Transmit spatial correlation matrix (or 3D array)
        %   Specify the spatial correlation of the transmitter as a double
        %   precision, real or complex, 2D matrix or 3D array. This
        %   property applies when you set the MIMOCorrelation property to
        %   'Custom' and the Polarization property to 'Co-Polar' or
        %   'Cross-Polar'. The first dimension of TransmitCorrelationMatrix
        %   should be the same as the number of transmit antennas Nt. If
        %   the channel is frequency-flat, i.e., PathDelays is a scalar,
        %   TransmitCorrelationMatrix is a 2D Hermitian matrix of size
        %   Nt-by-Nt. The main diagonal elements must be all ones, while
        %   the off-diagonal elements must be real or complex numbers with
        %   a magnitude smaller than or equal to one.
        %
        %   If the channel is frequency-selective, i.e., PathDelays is a
        %   row vector of length Np, TransmitCorrelationMatrix can be
        %   specified as a 2D matrix, in which case each path has the same
        %   transmit spatial correlation matrix. Alternatively, it can be
        %   specified as a 3-D array of size Nt-by-Nt-by-Np, in which case
        %   each path can have its own different transmit spatial
        %   correlation matrix.
        %
        %   The default value of this property is [1].
        TransmitCorrelationMatrix = 1;
        
        %ReceiveCorrelationMatrix Receive spatial correlation matrix (or 3D array)
        %   Specify the spatial correlation of the receiver as a double
        %   precision, real or complex, 2D matrix or 3D array. This
        %   property applies when you set the MIMOCorrelation property to
        %   'Custom' and the Polarization property to 'Co-Polar' or
        %   'Cross-Polar'. The first dimension of ReceiveCorrelationMatrix
        %   should be the same as the number of receive antennas Nr. If the
        %   channel is frequency-flat, i.e., PathDelays is a scalar,
        %   ReceiveCorrelationMatrix is a 2D Hermitian matrix of size
        %   Nr-by-Nr. The main diagonal elements must be all ones, while
        %   the off-diagonal elements must be real or complex numbers with
        %   a magnitude smaller than or equal to one.
        %  
        %   If the channel is frequency-selective, i.e., PathDelays is a
        %   row vector of length Np, ReceiveCorrelationMatrix can be
        %   specified as a 2D matrix, in which case each path has the same
        %   receive spatial correlation matrix. Alternatively, it can be
        %   specified as a 3-D array of size Nr-by-Nr-by-Np, in which case
        %   each path can have its own different receive spatial
        %   correlation matrix.
        % 
        %   The default value of this property is [1 0; 0 1].
        ReceiveCorrelationMatrix = eye(2);
        
        %TransmitPolarizationAngles Transmit polarization slant angles in degrees
        %   Specify the transmitter antenna polarization angles, in
        %   degrees, as a double-precision row vector. This property
        %   applies when MIMOCorrelation is set to 'Custom' and
        %   Polarization is set to 'Cross-Polar'.
        %
        %   The default value of this property is [45 -45].        
        TransmitPolarizationAngles = [45 -45];

        %ReceivePolarizationAngles Receive polarization slant angles in degrees
        %   Specify the receiver antenna polarization angles, in degrees,
        %   as a double-precision row vector. This property applies when
        %   MIMOCorrelation is set to 'Custom' and Polarization is set to
        %   'Cross-Polar'.
        %
        %   The default value of this property is [90 0].
        ReceivePolarizationAngles = [90 0];
        
        %XPR Cross polarization power ratio (dB)
        %   Specify the cross-polarization power ratio in dB as a scalar or
        %   row vector. The XPR is defined as used in the Clustered Delay
        %   Line (CDL) models in TR 38.901 Section 7.7.1, where the XPR is
        %   the ratio between the vertical-to-vertical and
        %   vertical-to-horizontal polarizations (P_vv/P_vh) i.e. the XPR
        %   in dB is zero or greater. This property applies when
        %   MIMOCorrelation is set to 'Custom' and Polarization is set to
        %   'Cross-Polar'.
        %
        %   If the channel is frequency-selective, i.e., PathDelays is a
        %   row vector of length Np, XPR can be specified as a scalar, in
        %   which case each path has the same XPR. Alternatively, it can be
        %   specified as a vector of size 1-by-Np, in which case each path
        %   can have its own different XPR.
        %
        %   The default value of this property is 10.0dB. 
        XPR = 10.0;

        %SpatialCorrelationMatrix Combined correlation matrix (or 3-D array)
        %   Specify the combined spatial correlation for the channel as a
        %   double precision, 2D matrix or 3D array. This property applies
        %   when you set the MIMOCorrelation property to 'Custom' and the
        %   Polarization property to 'Custom'. The first dimension of
        %   SpatialCorrelationMatrix determines the product of the number
        %   of transmit antennas Nt and the number of receive antennas Nr.
        %   If the channel is frequency-flat, i.e., PathDelays is a scalar,
        %   SpatialCorrelationMatrix is a 2D Hermitian matrix of size
        %   (NtNr)-by-(NtNr). The magnitude of any off-diagonal element
        %   must be no larger than the geometric mean of the two
        %   corresponding diagonal elements.
        %  
        %   If the channel is frequency-selective, i.e., PathDelays is a
        %   row vector of length Np, SpatialCorrelationMatrix can be
        %   specified as a 2D matrix, in which case each path has the same
        %   spatial correlation matrix. Alternatively, it can be specified
        %   as a 3-D array of size (NtNr)-by-(NtNr)-by-Np, in which case
        %   each path can have its own different spatial correlation
        %   matrix.
        % 
        %   The default value of this property is [1 0; 0 1].        
        SpatialCorrelationMatrix = eye(2);
        
    end
    
    properties (Access = public, Nontunable, Logical)
            
        %NormalizePathGains Normalize path gains to total power of 0 dB (logical)
        %   Set this property to true to normalize the fading processes
        %   such that the total power of the path gains, averaged over
        %   time, is 0 dB. When you set this property to false, there is no
        %   normalization on path gains. The average powers of the path
        %   gains are specified by the selected delay profile or by the
        %   AveragePathGains property if DelayProfile is set to 'Custom'.
        %
        %   The default value of this property is true. 
        NormalizePathGains = true;
                
        %PathGainsOutputPort Enable path gain output (logical)
        %   Set this property to true to output the channel path gains of
        %   the underlying fading process. 
        %
        %   The default value of this property is false.
        PathGainsOutputPort = false;
        
    end

    properties (Access = public, Nontunable)
        
        %InitialTime Start time of the fading process (s)
        %   Specify the time offset of the fading process as a real
        %   nonnegative scalar.
        %
        %   The default value of this property is 0.
        InitialTime = 0.0;
        
        %NumSinusoids Number of sinusoids used to model the fading process
        %   Specify the number of sinusoids used to model the channel as a
        %   positive integer scalar.
        %
        %   The default value of this property is 48.
        NumSinusoids = 48;
        
        %RandomStream Source of random number stream
        %   Specify the source of random number stream as one of 'Global
        %   stream' or 'mt19937ar with seed'. If you set RandomStream to
        %   'Global stream', the current global random number stream is
        %   used for normally distributed random number generation. In this
        %   case, the reset method only resets the filters. If you set
        %   RandomStream to 'mt19937ar with seed', the mt19937ar algorithm
        %   is used for normally distributed random number generation. In
        %   this case, the reset method not only resets the filters but
        %   also reinitializes the random number stream to the value of the
        %   Seed property.
        %
        %   The default value of this property is 'mt19937ar with seed'. 
        RandomStream = 'mt19937ar with seed';
        
        %Seed Initial seed of mt19937ar random number stream
        %   Specify the initial seed of a mt19937ar random number generator
        %   algorithm as a double precision, real, nonnegative integer
        %   scalar. This property applies when you set the RandomStream
        %   property to 'mt19937ar with seed'. The Seed reinitializes the
        %   mt19937ar random number stream in the reset method.
        %
        %   The default value of this property is 73.
        Seed = 73;
        
    end
    
    properties (Access = public, Nontunable, Logical)
        
        %NormalizeChannelOutputs Normalize channel outputs by the number of receive antennas (logical)
        %   Set this property to true to normalize the channel outputs by
        %   the number of receive antennas. When you set this property to
        %   false, there is no normalization for channel outputs.
        %
        %   The default value of this property is true.
        NormalizeChannelOutputs = true;
        
    end
    
    % public property setters for validation
    methods
        
        function set.PathDelays(obj,val)
            propName = 'PathDelays';
            validateattributes(val,{'double'},{'real','row','finite'},[class(obj) '.' propName],propName);
            obj.PathDelays = val;
        end
        
        function set.AveragePathGains(obj,val)
            propName = 'AveragePathGains';
            validateattributes(val,{'double'},{'real','row','finite'},[class(obj) '.' propName], propName);
            obj.AveragePathGains = val;
        end
        
        function set.KFactorFirstTap(obj,val)
            propName = 'KFactorFirstTap';
            validateattributes(val,{'double'},{'real','scalar','finite'},[class(obj) '.' propName],propName);
            obj.KFactorFirstTap = val;
        end
        
        function set.DelaySpread(obj,val)
            propName = 'DelaySpread';
            validateattributes(val,{'double'},{'real','scalar','nonnegative','finite'},[class(obj) '.' propName],propName);
            obj.DelaySpread = val;
        end
                
        function set.MaximumDopplerShift(obj,val)
            propName = 'MaximumDopplerShift';
            validateattributes(val,{'double'},{'real','scalar','nonnegative','finite'},[class(obj) '.' propName],propName);
            obj.MaximumDopplerShift = val;
        end
        
        function set.KFactor(obj,val)
            propName = 'KFactor';
            validateattributes(val,{'double'},{'real','scalar','finite'},[class(obj) '.' propName],propName);
            obj.KFactor = val;
        end
        
        function set.SampleRate(obj,val)
            propName = 'SampleRate';
            validateattributes(val,{'double'},{'real','scalar','positive','finite'},[class(obj) '.' propName],propName);
            obj.SampleRate = val;
        end
        
        function set.NumTransmitAntennas(obj,val)
            propName = 'NumTransmitAntennas';
            validateattributes(val,{'numeric'},{'real','scalar','integer','>=',1},[class(obj) '.' propName],propName);
            obj.NumTransmitAntennas = val;
        end
        
        function set.NumReceiveAntennas(obj,val)
            propName = 'NumReceiveAntennas';
            validateattributes(val,{'numeric'},{'real','scalar','integer','>=',1},[class(obj) '.' propName],propName);
            obj.NumReceiveAntennas = val;
        end
        
        function set.TransmitCorrelationMatrix(obj,val)
            propName = 'TransmitCorrelationMatrix';
            validateattributes(val,{'double'},{'finite','nonempty'},[class(obj) '.' propName],propName);
            coder.internal.errorIf(ndims(val)>3,'lte5g:nr5gTDLChannel:CorrMtxMoreThan3D','TransmitCorrelationMatrix');
            obj.TransmitCorrelationMatrix = val;
        end
        
        function set.ReceiveCorrelationMatrix(obj,val)
            propName = 'ReceiveCorrelationMatrix';
            validateattributes(val,{'double'},{'finite','nonempty'},[class(obj) '.' propName],propName);
            coder.internal.errorIf(ndims(val)>3,'lte5g:nr5gTDLChannel:CorrMtxMoreThan3D','ReceiveCorrelationMatrix');
            obj.ReceiveCorrelationMatrix = val;
        end
        
        function set.TransmitPolarizationAngles(obj,val)
            propName = 'TransmitPolarizationAngles';
            validateattributes(val,{'double'},{'real','row','finite'},[class(obj) '.' propName], propName);
            coder.internal.errorIf(size(val,2)>2,'lte5g:nr5gTDLChannel:InvalidNumPolAngles','TransmitPolarizationAngles');
            obj.TransmitPolarizationAngles = val;
        end
        
        function set.ReceivePolarizationAngles(obj,val)
            propName = 'ReceivePolarizationAngles';
            validateattributes(val,{'double'},{'real','row','finite'},[class(obj) '.' propName], propName);
            coder.internal.errorIf(size(val,2)>2,'lte5g:nr5gTDLChannel:InvalidNumPolAngles','ReceivePolarizationAngles');
            obj.ReceivePolarizationAngles = val;
        end
        
        function set.XPR(obj,val)
            propName = 'XPR';
            validateattributes(val,{'double'},{'real','row','nonnegative'},[class(obj) '.' propName],propName);
            obj.XPR = val;
        end
        
        function set.SpatialCorrelationMatrix(obj,val)
            propName = 'SpatialCorrelationMatrix';
            validateattributes(val,{'double'},{'finite','nonempty'},[class(obj) '.' propName],propName);
            coder.internal.errorIf(ndims(val)>3,'lte5g:nr5gTDLChannel:CorrMtxMoreThan3D','SpatialCorrelationMatrix');
            obj.SpatialCorrelationMatrix = val;
        end
        
        function set.InitialTime(obj,val)
            propName = 'InitialTime';
            validateattributes(val,{'double'},{'real','scalar','nonnegative','finite'},[class(obj) '.' propName],propName);
            obj.InitialTime = val;
        end
        
        function set.NumSinusoids(obj,val)
            propName = 'NumSinusoids';
            validateattributes(val,{'numeric'},{'scalar','integer','>=',1},[class(obj) '.' propName],propName);
            obj.NumSinusoids = val;
        end
        
        function set.Seed(obj,val)
            propName = 'Seed';
            validateattributes(val,{'double'},{'real','scalar','integer','nonnegative','finite'},[class(obj) '.' propName],propName);
            obj.Seed = val;
        end
        
    end
    
    % property value sets for enumerated properties
    properties(Hidden,Transient)
        
        DelayProfileSet = matlab.system.StringSet({'TDL-A','TDL-B','TDL-C','TDL-D','TDL-E','Custom'});
        MIMOCorrelationSet = matlab.system.StringSet({'Low','Medium','Medium-A','UplinkMedium','High','Custom'});
        FadingDistributionSet = matlab.system.StringSet({'Rayleigh','Rician'});
        RandomStreamSet = matlab.system.StringSet({'Global stream','mt19937ar with seed'});
        PolarizationSet = matlab.system.StringSet({'Co-Polar','Cross-Polar','Custom'});
        TransmissionDirectionSet = matlab.system.StringSet({'Downlink','Uplink'});
        
    end

% =========================================================================
%   protected interface
    
    methods (Access = protected)
        
        % nr5gTDLChannel setupImpl method
        function setupImpl(obj)
                        
            % Construct the MIMOChannel
            obj.theChannel = nr5gTDLChannel.makeMIMOChannel(obj);
            
        end
        
        % nr5gTDLChannel stepImpl method
        function varargout = stepImpl(obj,in)
            
            % execute the MIMOChannel
            [varargout{1:nargout}] = obj.theChannel(in);

        end
        
        % nr5gTDLChannel resetImpl method
        function resetImpl(obj)
            
            % reset the MIMOChannel
            reset(obj.theChannel);
            
        end
        
        % nr5gTDLChannel releaseImpl method
        function releaseImpl(obj)
            
            release(obj.theChannel);
            
        end
        
        % nr5gTDLChannel getNumOutputsImpl method
        function num = getNumOutputsImpl(obj)
            
            num = 1 + obj.PathGainsOutputPort;
            
        end
        
        % nr5gTDLChannel infoImpl method
        function s = infoImpl(obj)
        %info Returns characteristic information about the TDL channel
        %   S = info(OBJ) returns a structure containing characteristic
        %   information, S, about the TDL fading channel. A description of
        %   the fields and their values is as follows:
        % 
        %   ChannelFilterDelay       - Channel filter delay in samples.
        %   AveragePathGains         - A row vector of the average gains of the
        %                              discrete paths, in dB. These values
        %                              include the effect of K-factor scaling if
        %                              enabled. 
        %   PathDelays               - A row vector providing the delays of the
        %                              discrete channel paths, in seconds. These
        %                              values include the effect of the desired
        %                              delay spread scaling, and desired
        %                              K-factor scaling if enabled.
        %   KFactorFirstTap          - K-factor of first tap of delay profile,
        %                              in dB. If the first tap of the delay
        %                              profile follows a Rayleigh rather than
        %                              Rician distribution, KFactorFirstTap will
        %                              be -Inf.
        %   NumTransmitAntennas      - Number of transmit antennas.
        %   NumReceiveAntennas       - Number of receive antennas.
        %   SpatialCorrelationMatrix - Combined correlation matrix (or 3-D array).
            
            if (~isempty(coder.target) || ~isLocked(obj))
                nr5gTDLChannel.validation(obj);
                channel = nr5gTDLChannel.makeMIMOChannel(obj);
            else
                channel = obj.theChannel;
            end
            
            theInfo = info(channel);
            s.ChannelFilterDelay = theInfo.ChannelFilterDelay;
            s.PathDelays = channel.PathDelays;
            s.AveragePathGains = channel.AveragePathGains;
            if (nr5gTDLChannel.hasLOSPath(obj))
                s.KFactorFirstTap = 10*log10(channel.KFactor);
            else
                s.KFactorFirstTap = -Inf;
            end
            Nt = channel.NumTransmitAntennas;
            s.NumTransmitAntennas = Nt;
            s.NumReceiveAntennas = size(channel.SpatialCorrelationMatrix,1) / Nt;
            s.SpatialCorrelationMatrix = channel.SpatialCorrelationMatrix;
            
        end
        
        % nr5gTDLChannel saveObjectImpl method
        function s = saveObjectImpl(obj)
            
            s = saveObjectImpl@matlab.System(obj);
            s.theChannel = matlab.System.saveObject(obj.theChannel);
            
        end

        % nr5gTDLChannel loadObjectImpl method
        function loadObjectImpl(obj,s,wasLocked)
            
            obj.theChannel = matlab.System.loadObject(s.theChannel);
            loadObjectImpl@matlab.System(obj,s,wasLocked);
            
        end
        
        % nr5gTDLChannel isInactivePropertyImpl method
        function flag = isInactivePropertyImpl(obj,prop)
            
            if (any(strcmp(prop,{'PathDelays','AveragePathGains','FadingDistribution'})))
                flag = ~strcmp(obj.DelayProfile,'Custom');
            elseif (strcmp(prop,'KFactorFirstTap'))
                flag = ~strcmp(obj.DelayProfile,'Custom') || ~strcmp(obj.FadingDistribution,'Rician');
            elseif (strcmp(prop,'KFactorScaling'))
                flag = ~any(strcmp(obj.DelayProfile,{'TDL-D','TDL-E'}));
            elseif (any(strcmp(prop,{'TransmitCorrelationMatrix','ReceiveCorrelationMatrix'})))
                flag = ~strcmp(obj.MIMOCorrelation,'Custom') || strcmp(obj.Polarization,'Custom');
            elseif (strcmp(prop,'NumTransmitAntennas'))
                flag = strcmp(obj.MIMOCorrelation,'Custom') && ~strcmp(obj.Polarization,'Custom');
            elseif (any(strcmp(prop,{'NumReceiveAntennas','TransmissionDirection'})))
                flag = strcmp(obj.MIMOCorrelation,'Custom');
            elseif (any(strcmp(prop,{'TransmitPolarizationAngles','ReceivePolarizationAngles','XPR'})))
                flag = ~(strcmp(obj.MIMOCorrelation,'Custom') && strcmp(obj.Polarization,'Cross-Polar'));
            elseif (strcmp(prop,'SpatialCorrelationMatrix'))
                flag = ~(strcmp(obj.MIMOCorrelation,'Custom') && strcmp(obj.Polarization,'Custom'));
            elseif (strcmp(prop,'KFactor'))
                flag = strcmp(obj.DelayProfile,'Custom') || ~(obj.KFactorScaling);
            elseif (strcmp(prop,'Seed'))
                flag = ~strcmp(obj.RandomStream,'mt19937ar with seed');
            elseif (strcmp(prop,'DelaySpread'))
                flag = strcmp(obj.DelayProfile,'Custom');
            else
                flag = false;
            end
            
        end
        
        % nr5gTDLChannel validatePropertiesImpl method
        function validatePropertiesImpl(obj)
            
            nr5gTDLChannel.validation(obj);
            
        end
    
    end

% =========================================================================
%   private

    properties (Access = private, Nontunable)
        
        % the underlying System object used to perform the channel modeling
        theChannel;
        
    end
    
    methods (Static, Access = private)
        
        function pathFilters = makePathFilters(obj)
            
            if (~isempty(coder.target) || ~isLocked(obj))
                nr5gTDLChannel.validation(obj);
                channel = nr5gTDLChannel.makeMIMOChannel(obj);
            else
                channel = obj.theChannel;
            end
            
            channelInfo = info(channel);
            
            pathFilters = channelInfo.ChannelFilterCoefficients.';
            
        end
        
        function c = makeMIMOChannel(obj)
            
            c = comm.MIMOChannel;
            c.SampleRate = obj.SampleRate;
            [pathDelays,averagePathGains,K_1dB] = nr5gTDLChannel.getDelayProfile(obj);
            c.PathDelays = pathDelays;
            c.AveragePathGains = averagePathGains;
            c.NormalizePathGains = obj.NormalizePathGains;
            if (nr5gTDLChannel.hasLOSPath(obj))
                c.FadingDistribution = 'Rician';
                c.KFactor = 10^(K_1dB/10); % convert dB to linear
                c.DirectPathDopplerShift = 0.7*obj.MaximumDopplerShift;
                c.DirectPathInitialPhase = 0.0;
            else
                c.FadingDistribution = 'Rayleigh';
            end
            c.MaximumDopplerShift = obj.MaximumDopplerShift;
            c.DopplerSpectrum = doppler('Jakes');
            c.SpatialCorrelationSpecification = 'Combined';
            c.NumTransmitAntennas = nr5gTDLChannel.getAntennaCount(obj);
            % NumReceiveAntennas is not used as SpatialCorrelationSpecification = 'Combined'
            % TransmitCorrelationMatrix is not used as SpatialCorrelationSpecification = 'Combined'
            % ReceiveCorrelationMatrix is not used as SpatialCorrelationSpecification = 'Combined'
            c.SpatialCorrelationMatrix = nr5gTDLChannel.getSpatialCorrelationMatrix(obj);
            c.AntennaSelection = 'Off';
            c.NormalizeChannelOutputs = obj.NormalizeChannelOutputs;
            c.FadingTechnique = 'Sum of sinusoids';
            c.NumSinusoids = obj.NumSinusoids;
            c.InitialTimeSource = 'Property';
            c.InitialTime = obj.InitialTime;
            c.RandomStream = obj.RandomStream;
            c.Seed = obj.Seed;
            c.PathGainsOutputPort = obj.PathGainsOutputPort;
            c.Visualization = 'Off';
            
        end
        
        function Rspat = getSpatialCorrelationMatrix(obj)
            
            coder.extrinsic('nr5g.internal.calculateSpatialCorrelationMatrix')
            
            if (strcmp(obj.MIMOCorrelation,'Custom') && strcmp(obj.Polarization,'Custom'))
            
                % Use SpatialCorrelationMatrix property as the overall
                % spatial correlation matrix Rspat
                Rspat = obj.SpatialCorrelationMatrix;
                
            else
                
                % Note: where the text "(-by-Np)" appears after a matrix
                % size here, this means that the matrix may instead be a
                % 3-D array with 3rd dimension size Np, the number of paths
                % i.e. the length of the PathDelays property (each path can
                % have its own different matrix)
                
                % Get permutation matrix according to TS 36.101 Annex
                % B.2.3A.1 / TS 36.104 Annex B.5A.1, size
                % (Nt*Nr)-by-(Nt*Nr) for cross-polar antennas. For co-polar
                % antennas, P is an identity matrix and therefore no
                % permutation is performed
                P = nr5gTDLChannel.getPermutationMatrix(obj);
                
                % Get transmit and receive side spatial correlation
                % matrices. For co-polar antennas, Rt and Rr are of size
                % Nt-by-Nt(-by-Np) and Nr-by-Nr(-by-Np) respectively. For
                % cross-polar antennas, they are of size
                % (Nt/2)-by-(Nt/2)(-by-Np) and (Nr/2)-by-(Nr/2)(-by-Np)
                % unless Nt=1 or Nr=1 in which case the size is
                % 1-by-1(-by-Np)
                Rt = nr5gTDLChannel.getTransmitCorrelationMatrix(obj);
                Rr = nr5gTDLChannel.getReceiveCorrelationMatrix(obj);
                
                % Get polarization covariance matrix. For co-polar antennas
                % Gamma is unity (size 1-by-1) and for cross-polar antennas
                % it is of size 4-by-4(-by-Np) unless Nt=1 or Nr=1 in which
                % case it is of size 2-by-2(-by-Np). (If both Nt and Nr are
                % 1, Gamma is of size 1-by-1(-by-Np))
                Gamma = nr5gTDLChannel.getPolarizationCorrelationMatrix(obj);
                
                % Compute overall spatial correlation matrix. Rspat is of
                % size (Nt*Nr)-by-(Nt*Nr)(-by-Np)
                a = nr5gTDLChannel.getRoundoffScalingFactor(obj);
                Rspat = coder.const(double(nr5g.internal.calculateSpatialCorrelationMatrix(P,Rt,Rr,Gamma,a)));       
                
            end
        
        end
        
        % TS 36.101 Annex B.2.3A.1
        % TS 36.104 Annex B.5A.1
        function P = getPermutationMatrix(obj)
            
            [Nt,Nr] = nr5gTDLChannel.getAntennaCount(obj);
            
            if (strcmp(obj.Polarization,'Co-Polar'))
                
                P = eye(Nr*Nt);
                
            else
                
                P = zeros(Nr*Nt);

                for i = 1:Nr
                    for j = 1:Nt

                        a = ((j-1)*Nr) + i;

                        if (j <= Nt/2)
                            b = (2*(j-1)*Nr) + i;
                        else
                            b = (2*(j-(Nt/2))*Nr) - Nr + i;
                        end

                        P(a,b) = 1;

                    end
                end
                
            end
            
        end
        
        function [Nt,Nr] = getAntennaCount(obj)
            
            if (strcmp(obj.MIMOCorrelation,'Custom'))
                if (strcmp(obj.Polarization,'Custom'))
                    Nt = obj.NumTransmitAntennas;
                    Nr = size(obj.SpatialCorrelationMatrix,1) / Nt;
                else
                    Nt_in = size(obj.TransmitCorrelationMatrix,1);
                    Nr_in = size(obj.ReceiveCorrelationMatrix,1);
                    if (strcmp(obj.Polarization,'Cross-Polar'))
                        [txPolAngles,rxPolAngles] = nr5gTDLChannel.getPolarizationAngles(obj);
                        Npt = numel(txPolAngles);
                        Npr = numel(rxPolAngles);
                        Nt = Nt_in * Npt;
                        Nr = Nr_in * Npr;
                    else
                        Nt = Nt_in;
                        Nr = Nr_in;
                    end
                end
            else
                Nt = obj.NumTransmitAntennas;
                Nr = obj.NumReceiveAntennas;
            end
            
        end
        
        % TS 36.101 Table B.2.3.2-1
        % TS 36.101 Table B.2.3A.3-1
        % TS 36.104 Table B.5.2-1
        % TS 36.104 Table B.5A.3-1
        function Rt = getTransmitCorrelationMatrix(obj)
            
            if (strcmpi(obj.MIMOCorrelation,'Custom'))
                Rt = obj.TransmitCorrelationMatrix;
            else
                switch (nr5gTDLChannel.getMIMOCorrelationVersusLink(obj))
                    case 'Low'
                        r = 0;
                    case {'Medium','Medium-A','UplinkMedium'}
                        r = 0.3;
                    case 'High'
                        r = 0.9;
                end
                Rt = nr5gTDLChannel.getCorrelationMatrix(obj,obj.NumTransmitAntennas,r);
            end

        end
        
        % TS 36.101 Table B.2.3.2-1
        % TS 36.101 Table B.2.3A.3-1
        % TS 36.104 Table B.5.2-1
        % TS 36.104 Table B.5A.3-1
        function Rr = getReceiveCorrelationMatrix(obj)

            if (strcmpi(obj.MIMOCorrelation,'Custom'))
                Rr = obj.ReceiveCorrelationMatrix;
            else
                switch (nr5gTDLChannel.getMIMOCorrelationVersusLink(obj))
                    case 'Low'
                        r = 0;
                    case {'Medium','UplinkMedium'}
                        r = 0.9;
                    case 'Medium-A'
                        if (strcmpi(obj.Polarization,'Co-Polar'))
                            r = 0.3874;
                        else
                            r = 0.6;
                        end
                    case 'High'
                        r = 0.9;
                end
                Rr = nr5gTDLChannel.getCorrelationMatrix(obj,obj.NumReceiveAntennas,r);
            end

        end
        
        function mimoCorrelation = getMIMOCorrelationVersusLink(obj)
            
            if (strcmpi(obj.TransmissionDirection,'Uplink') && strcmpi(obj.MIMOCorrelation,'Medium'))
                mimoCorrelation = 'UplinkMedium';
            else
                mimoCorrelation = obj.MIMOCorrelation;
            end
            
        end
                
        function Gamma = getPolarizationCorrelationMatrix(obj)
            
            coder.extrinsic('nr5g.internal.calculateGammaMatrix');
            
            switch (obj.Polarization)
                
                case 'Co-Polar'
                    
                    % For co-polar antennas, correlation is defined only in
                    % terms of spatial positions so polarization
                    % correlation matrix is unity (size 1-by-1)
                    Gamma = 1;
                    
                case 'Cross-Polar'
                    
                    % Transmit and receive polarization angles
                    [txPolAngles,rxPolAngles] = nr5gTDLChannel.getPolarizationAngles(obj);
                    
                    % Cross polarization power ratio in dB
                    if (strcmp(obj.MIMOCorrelation,'Custom'))
                        XPR = -obj.XPR;
                    else
                        gamma = nr5gTDLChannel.getPolarizationCorrelation(obj);
                        XPR = 10*log10((1-gamma)/(1+gamma));
                    end
                    
                    % Calculate polarization correlation matrix according
                    % to IEEE 802.16m-08/004r5 Appendix B
                    Gamma = coder.const(double(nr5g.internal.calculateGammaMatrix(txPolAngles,rxPolAngles,XPR)));
                    
            end
            
        end
        
        function [txPolAnglesOut,rxPolAnglesOut] = getPolarizationAngles(obj)
            
            if (strcmp(obj.MIMOCorrelation,'Custom'))
                txPolAnglesOut = obj.TransmitPolarizationAngles;
                rxPolAnglesOut = obj.ReceivePolarizationAngles;
            else
                if (strcmp(obj.Polarization,'Co-Polar'))
                    % TS 36.101 Annex B.2.3
                    % TS 36.104 Annex B.5
                    txPolAnglesOut = 90;
                    rxPolAnglesOut = 90;
                else
                    % TS 36.101 Annex B.2.3A
                    % TS 36.104 Annex B.5A
                    if (strcmp(obj.TransmissionDirection,'Downlink'))
                        txPolAngles = [45 -45];
                        rxPolAngles = [90 0];
                    else
                        txPolAngles = [90 0];
                        rxPolAngles = [45 -45];
                    end
                    if (obj.NumTransmitAntennas==1)
                        txPolAnglesOut = txPolAngles(1);
                    else
                        txPolAnglesOut = txPolAngles;
                    end
                    if (obj.NumReceiveAntennas==1)
                        rxPolAnglesOut = rxPolAngles(1);
                    else
                        rxPolAnglesOut = rxPolAngles;
                    end
                end
            end

        end
            
        % TS 36.101 Table B.2.3A.3-1
        % TS 36.104 Table B.5A.3-1
        function gamma = getPolarizationCorrelation(obj)
           
            switch (nr5gTDLChannel.getMIMOCorrelationVersusLink(obj))
                case 'Low'
                    gamma = 0;
                case {'Medium','Medium-A','UplinkMedium'}
                    gamma = 0.2;
                case 'High'
                    gamma = 0.3;
            end
            
        end

        % TS 36.101 Annex B.2.3.2
        % TS 36.101 Annex B.2.3A.3
        % TS 36.104 Annex B.5.2
        function a = getRoundoffScalingFactor(obj)
            
            [Nt,Nr] = nr5gTDLChannel.getAntennaCount(obj);
            
            a = 0.0;
            if (~strcmpi(obj.MIMOCorrelation,'Custom'))
                if (Nt==4 && Nr==4)
                    a = 0.00012;
                elseif (strcmpi(obj.MIMOCorrelation,'High') && ((((Nt==2 && Nr==4) || (Nt==4 && Nr==2)) && strcmpi(obj.Polarization,'Co-Polar')) || (Nt==8 && Nr==2)))
                    a = 0.00010;
                end
            end
            
        end
        
        % TS 36.101 Annex B.2.3.1
        % TS 36.101 Annex B.2.3A.2
        % TS 36.104 Annex B.5.1
        % TS 36.104 Annex B.5A.2
        function R = getCorrelationMatrix(obj,N,r)

            if (strcmp(obj.Polarization,'Cross-Polar') && N~=1)
                N = N/2;
            end
                
            switch (N)
                case 1
                    R = 1;
                case 2
                    R = [1 r; ...
                         r 1];
                case 4
                    R = [1          r^(1/9)    r^(4/9)    r;       ...
                         (r^(1/9))  1          r^(1/9)    r^(4/9); ...
                         (r^(4/9))  (r^(1/9))  1          r^(1/9); ...
                         r          (r^(4/9))  (r^(1/9))  1];
                otherwise
                    R = eye(N);
            end

        end

        function [pathDelaysOut,pathGainsOut,K_1dB] = getDelayProfile(obj)
            
            coder.extrinsic('nr5g.internal.getTDLProfile');
            coder.extrinsic('lte.internal.scaleDelaysAndKFactor');
            
            if(strcmpi(obj.DelayProfile,'Custom'))
                pathDelaysIn = obj.PathDelays;
                pathGainsIn = obj.AveragePathGains;
                if (nr5gTDLChannel.hasLOSPath(obj))
                    % Split the first path into a LOS part and Rayleigh
                    % part according to K_1
                    K_1dB = obj.KFactorFirstTap;
                    K_1 = 10^(K_1dB/10);
                    P_1dB = pathGainsIn(1);
                    P_1 = 10^(P_1dB/10);
                    pathDelays = [pathDelaysIn(1) pathDelaysIn(1) pathDelaysIn(2:end)];
                    pathGains = [(-10*log10(P_1 + P_1/K_1) + [0 -K_1dB]) pathGainsIn(2:end)];
                else
                    pathDelays = pathDelaysIn;
                    pathGains = pathGainsIn;
                end
            else
                desiredKFactor = NaN;
                if (nr5gTDLChannel.hasLOSPath(obj))
                    if (obj.KFactorScaling)
                        desiredKFactor = obj.KFactor;
                    end
                end
                pdp = coder.const(double(nr5g.internal.getTDLProfile(obj.DelayProfile)));
                pdp = coder.const(double(lte.internal.scaleDelaysAndKFactor(pdp,desiredKFactor,obj.DelaySpread)));
                pathDelays = pdp(:,1).'; % 1st column is delay
                pathGains = pdp(:,2).';  % 2nd column is power
            end
            
            % At this point in the code, if a Rician path is present it is
            % split into a LOS part and a Rayleigh part, whether the delay
            % profile was from the standard tables or is custom
            
            if (nr5gTDLChannel.hasLOSPath(obj))
                % Remove 2nd entry of delay profile, corresponding to the
                % Rayleigh part of the Rician path. The first entry plus
                % K_1dB now capture the Rician path
                K_1dB = pathGains(1) - pathGains(2);
                pathGainsOut = [10*log10(sum(10.^(pathGains(1:2)/10))) pathGains(3:end)];
                pathDelaysOut = pathDelays([1 3:end]);
            else
                K_1dB = -Inf;
                pathGainsOut = pathGains;
                pathDelaysOut = pathDelays;
            end
            
        end

        function has = hasLOSPath(obj)

            if (strcmpi(obj.DelayProfile,'Custom'))
                has = strcmpi(obj.FadingDistribution,'Rician');
            else
                has = lte.internal.hasLOSPath(obj.DelayProfile);
            end

        end
        
        function validation(obj)
            
            coder.internal.errorIf( ...
                strcmp(obj.Polarization,'Custom') && ~strcmp(obj.MIMOCorrelation,'Custom'), ...
                'lte5g:nr5gTDLChannel:InvalidCorrelationForPolarization');
            
            [txPolAngles,rxPolAngles] = nr5gTDLChannel.getPolarizationAngles(obj);
            Npt = numel(txPolAngles);
            Npr = numel(rxPolAngles);
            
            coder.internal.errorIf( ...
                ~strcmp(obj.MIMOCorrelation,'Custom') && ...
                strcmp(obj.Polarization,'Cross-Polar') && ...
                (mod(obj.NumTransmitAntennas,Npt)~=0 || mod(obj.NumReceiveAntennas,Npr)~=0), ...
                'lte5g:nr5gTDLChannel:InvalidAntsForCrossPolar');
            
            coder.internal.errorIf( ...
                any(strcmp(obj.MIMOCorrelation, {'Medium', 'Medium-A', 'High'})) && ...
                (all(obj.NumTransmitAntennas ~= [1 2 4]*Npt) || ...
                 all(obj.NumReceiveAntennas  ~= [1 2 4]*Npr)), ...
                'lte5g:nr5gTDLChannel:InvalidAntsForCorrelation');
            
            Np = size(nr5gTDLChannel.getDelayProfile(obj),2);
            nr5gTDLChannel.validateCorrelationMatrixSize(nr5gTDLChannel.getTransmitCorrelationMatrix(obj),Np,'TransmitCorrelationMatrix');
            nr5gTDLChannel.validateCorrelationMatrixSize(nr5gTDLChannel.getReceiveCorrelationMatrix(obj),Np,'ReceiveCorrelationMatrix');
            if (strcmp(obj.MIMOCorrelation,'Custom') && strcmp(obj.Polarization,'Cross-Polar'))
                coder.internal.errorIf(all(size(obj.XPR,2)~=[1 Np]),'lte5g:nr5gTDLChannel:XPRDimNotMatchNP');
            end
            
        end
        
        function validateCorrelationMatrixSize(R,Np,name)
            
            coder.internal.errorIf(ndims(R)==3 && size(R,3)~=Np,'lte5g:nr5gTDLChannel:CorrMtxDimNotMatchNP',name);
            
        end

    end
    
end
