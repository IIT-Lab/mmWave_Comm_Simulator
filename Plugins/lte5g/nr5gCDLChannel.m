classdef (StrictDefaults)nr5gCDLChannel < matlab.System
%nr5gCDLChannel TR 38.901 Clustered Delay Line (CDL) channel
%   CHAN = nr5gCDLChannel creates a CDL MIMO fading channel System object,
%   CHAN. This object filters a real or complex input signal through the
%   CDL MIMO channel to obtain the channel impaired signal. This object
%   implements the following aspects of TR 38.901:
%   * Section 7.7.1 Clustered Delay Line (CDL) models
%   * Section 7.7.3 Scaling of delays 
%   * Section 7.7.6 K-factor for LOS channel models
%   * Section 7.7.5.1 Scaling of angles
%
%   CHAN = nr5gCDLChannel(Name,Value) creates a CDL MIMO channel object,
%   CHAN, with the specified property Name set to the specified Value. You
%   can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax for ChannelFiltering = true:
%
%   Y = step(CHAN,X) filters input signal X through a CDL MIMO fading
%   channel and returns the result in Y. The input X can be a double or
%   single precision data type scalar, vector, or 2D matrix with real or
%   complex values. X is of size Ns-by-Nt, where Ns is the number of
%   samples and Nt is the number of transmit antennas. Y is the output
%   signal of size Ns-by-Nr, where Nr is the number of receive antennas. Y
%   is of double or single precision data type with complex values, where
%   the precision is the same as for input X. 
% 
%   [Y,PATHGAINS] = step(CHAN,X) also returns the MIMO channel path gains
%   of the underlying fading process in PATHGAINS. PATHGAINS is of size
%   Ncs-by-Np-by-Nt-by-Nr, where Np is the number of paths, i.e., the
%   length of the PathDelays property value of CHAN, and Ncs is the number
%   of channel snapshots, controlled by the SampleDensity property.
%   PATHGAINS is of the same precision as input X with complex values.
%
%   [Y,PATHGAINS,SAMPLETIMES] = step(CHAN,X) also returns the sample times
%   of the channel snapshots (1st dimension elements) of PATHGAINS.
%   SAMPLETIMES is of size Ncs-by-1 and is of double precision data type
%   with real values.
%
%   Step method syntax for ChannelFiltering = false:
%
%   [PATHGAINS,SAMPLETIMES] = step(CHAN) produces path gains PATHGAINS and
%   sample times SAMPLETIMES as described above, where the duration of the
%   fading process is given by the NumTimeSamples property. In this case
%   the object acts as a source of path gains and sample times without
%   filtering an input signal. 
% 
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(obj,x) and y = obj(x) are
%   equivalent.
%
%   nr5gCDLChannel methods:
%
%   step           - Filter input signal through a CDL MIMO fading channel
%                    (see above)
%   release        - Allow property value and input characteristics changes
%   clone          - Create a CDL channel object with same property values
%   isLocked       - Locked status (logical)
%   reset          - Reset states of filters, and random stream if the
%                    RandomStream property is set to 'mt19937ar with seed'
%   info           - Return characteristic information about the CDL 
%                    channel
%   getPathFilters - Get filter impulse responses for the filters which 
%                    apply the path delays to the input waveform
%
%   nr5gCDLChannel properties:
%
%   DelayProfile            - CDL delay profile
%   PathDelays              - Discrete path delay vector (s)
%   AveragePathGains        - Average path gain vector (dB)
%   AnglesAoD               - Azimuth of departure angles vector (deg)
%   AnglesAoA               - Azimuth of arrival angles vector (deg)
%   AnglesZoD               - Zenith of departure angles vector (deg)
%   AnglesZoA               - Zenith of arrival angles vector (deg)
%   HasLOSCluster           - Line of sight cluster (logical)
%   KFactorFirstCluster     - K-factor of first cluster (dB)
%   AngleScaling            - Enable angle scaling (logical)
%   AngleSpreads            - Desired scaled angle spreads vector (deg)
%   MeanAngles              - Desired mean scaled angles vector (deg)
%   XPR                     - Cross polarization power ratio (dB)
%   DelaySpread             - Desired delay spread (s)
%   CarrierFrequency        - Carrier frequency (Hz)
%   MaximumDopplerShift     - Maximum Doppler shift (Hz)
%   UTDirectionOfTravel     - User terminal (UT) direction of travel
%   KFactorScaling          - Enable K-factor scaling (logical)
%   KFactor                 - Desired Rician K-factor (dB)
%   SampleRate              - Input signal sample rate (Hz)
%   TransmitAntennaArray    - Transmit antenna array characteristics
%   ReceiveAntennaArray     - Receive antenna array characteristics
%   SampleDensity           - Number of time samples per half wavelength 
%   NormalizePathGains      - Normalize path gains (logical)
%   InitialTime             - Start time of fading process (s)
%   NumStrongestClusters    - Number of strongest clusters to split into subclusters
%   ClusterDelaySpread      - Cluster delay spread (s)
%   RandomStream            - Source of random number stream
%   Seed                    - Initial seed of mt19937ar random number stream
%   ChannelFiltering        - Perform filtering of input signal (logical)
%   NumTimeSamples          - Number of time samples
%   NormalizeChannelOutputs - Normalize channel outputs (logical)
% 
%   % Example 1:
%   %   Create an LTE waveform for Reference Measurement Channel (RMC) R.50
%   %   TDD (10MHz, QPSK, R=1/3, 1 layer, 8 CSI-RS ports) and pass it
%   %   through a channel with delay profile CDL-D, 10ns delay spread and
%   %   UT velocity 15km/h:
%
%   rmc = lteRMCDL('R.50','TDD');
%   rmc.TotSubframes = 1;
%   data = [1; 0; 0; 1];
%   [txWaveform,~,txInfo] = lteRMCDLTool(rmc,data);
%
%   v = 15.0;                    % UT velocity in km/h
%   fc = 4e9;                    % carrier frequency in Hz
%   c = physconst('lightspeed'); % speed of light in m/s
%   fd = (v*1000/3600)/c*fc;     % UT max Doppler frequency in Hz
%
%   cdl = nr5gCDLChannel;
%   cdl.DelayProfile = 'CDL-D';
%   cdl.DelaySpread = 10e-9;
%   cdl.CarrierFrequency = fc;
%   cdl.MaximumDopplerShift = fd;
%   cdl.SampleRate = txInfo.SamplingRate;
%
%   % Configure transmit and receive antenna arrays. The transmit array is
%   % configured as [M N P M_g N_g] = [2 2 2 1 1], i.e. 1 panel (M_g=1,
%   % N_g=1) with 2x2 antennas (M=2, N=2) and P=2 polarization angles. The
%   % receive antenna array is configured as [M N P M_g N_g] = [1 1 2 1 1]
%   % i.e. a single pair of cross-polarized co-located antennas
%
%   cdl.TransmitAntennaArray.Size = [2 2 2 1 1];
%   cdl.ReceiveAntennaArray.Size = [1 1 2 1 1];
%
%   rxWaveform = cdl(txWaveform);
%
%   % Example 2: 
%   %   Configure a channel for SISO operation and delay profile CDL-B. Set
%   %   the maximum Doppler shift to 300Hz and the channel sampling rate to
%   %   10kHz. Then plot the step response of the channel (lines) and the 
%   %   corresponding path gain snapshots (circles) for various values of 
%   %   the SampleDensity property, which controls how often the channel 
%   %   snapshots are taken relative to the Doppler frequency. 
%   %   SampleDensity=Inf ensures that a channel snapshot is taken for 
%   %   every input sample. SampleDensity=X takes channel snapshots at a 
%   %   rate of Fcs=2*X*MaximumDopplerShift. The channel snapshots are 
%   %   applied to the input waveform by means of zero order hold 
%   %   interpolation. Note that an extra snapshot is taken beyond the end 
%   %   of the input, where some of the final output samples use this extra 
%   %   value to help minimize the interpolation error. Note that the 
%   %   channel output contains a transient (and delay) due to the filters
%   %   that implement the path delays.
%
%   cdl = nr5gCDLChannel;
%   cdl.TransmitAntennaArray.Size = [1 1 1 1 1];
%   cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
%   cdl.DelayProfile = 'CDL-B';
%   cdl.MaximumDopplerShift = 300.0;
%   cdl.SampleRate = 10e3;
%   cdl.Seed = 19;
%   T = 40; in = ones(T,1); SR = cdl.SampleRate;
%   disp(['input length T=' num2str(T) ' samples']);
% 
%   s = [Inf 5 2]; % sample densities
% 
%   legends = {};
%   figure; hold on;
%   for i = 1:length(s)
%     
%       % execute channel with chosen sample density
%       release(cdl); cdl.SampleDensity = s(i);
%       [out,pathgains,sampletimes] = cdl(in);
%       chInfo = info(cdl); tau = chInfo.ChannelFilterDelay;
%     
%       % plot channel output against time
%       t = cdl.InitialTime + ((0:(T-1)) - tau).' / SR;
%       h = plot(t,abs(out),'o-'); h.MarkerSize = 2; h.LineWidth = 1.5;
%       desc = ['SampleDensity=' num2str(s(i))];
%       legends = [legends ['output, ' desc]];
%       disp([desc ', Ncs=' num2str(length(sampletimes))]);
%     
%       % plot path gains against sample times
%       h2 = plot(sampletimes-tau/SR,abs(sum(pathgains,2)),'o');
%       h2.Color = h.Color; h2.MarkerFaceColor = h.Color;
%       legends = [legends ['path gains, ' desc]];
%     
%   end
%   xlabel('time (s)');
%   title('Channel output and path gains versus SampleDensity');
%   ylabel('channel magnitude');
%   legend(legends,'Location','NorthWest');
%
%   % Example 3:
%   %   Create a resource grid for 64 antennas, fill it with QPSK symbols
%   %   and perform LTE OFDM modulation. Pass the resulting waveform
%   %   through a 64-by-4 channel and plot the received waveform spectrum:
%
%   enb.NDLRB = 25;
%   enb.CyclicPrefix = 'Normal';
%   grid = lteDLResourceGrid(enb,64);
%   grid(:) = lteSymbolModulate(randi([0 1],numel(grid)*2,1),'QPSK');
%   [txWaveform,txInfo] = lteOFDMModulate(enb,grid);
%
%   cdl = nr5gCDLChannel;
%   cdl.SampleRate = txInfo.SamplingRate;
%   cdl.TransmitAntennaArray.Size = [2 4 2 2 2];
%   cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 2.0 1.0];
%   cdl.ReceiveAntennaArray.Size = [2 1 2 1 1];
%   
%   % The antenna array elements are mapped to the waveform channels
%   % (columns) in the order that a 5-D array of size
%   % TransmitAntennaArray.Size or ReceiveAntennaArray.Size is linearly
%   % indexed i.e. across the dimensions first to last. See the
%   % TransmitAntennaArray or ReceiveAntennaArray property help for more
%   % details.
%   
%   rxWaveform = cdl(txWaveform);
%
%   analyzer = dsp.SpectrumAnalyzer('SampleRate',cdl.SampleRate);
%   analyzer.Title = ['Received signal spectrum for ' cdl.DelayProfile];
%   analyzer(rxWaveform);
%
%   See also nr5gTDLChannel, lteFadingChannel, comm.MIMOChannel.

%   Copyright 2016-2018 The MathWorks, Inc.
    
%#codegen

% =========================================================================
%   public interface

    methods (Access = public)
        
        % nr5gCDLChannel constructor
        function obj = nr5gCDLChannel(varargin)
            
            % Set property values from any name-value pairs input to the
            % constructor
            setProperties(obj,nargin,varargin{:});
            
        end
        
        function h = getPathFilters(obj)
        %getPathFilters Get path filter impulse responses
        %   H = getPathFilters(obj) returns a double precision real matrix
        %   of size Nh-by-Np where Nh is the number of impulse response
        %   samples and Np is the number of paths. Each column of H
        %   contains the filter impulse response for each path of the delay
        %   profile. This information facilitates reconstruction of a
        %   perfect channel estimate when used in conjunction with the
        %   PATHGAINS output of the step method.
            
            if (isempty(coder.target) && obj.ChannelFiltering && isLocked(obj))
                h = obj.pathFilters;
            else
                objinfo = info(obj);
                h = lte.internal.designCDLPathDelayFilters(objinfo.PathDelays,obj.channelFilterDelay,obj.maxFractionalDelayError,obj.SampleRate);
            end
            
        end
        
    end
    
    properties (Access = public, Nontunable)
        
        %DelayProfile CDL delay profile
        %   Specify the CDL delay profile as one of 'CDL-A', 'CDL-B',
        %   'CDL-C', 'CDL-D', 'CDL-E' or 'Custom'. See TR 38.901 Section
        %   7.7.1, Tables 7.7.1-1 to 7.7.1-5. When you set this property to
        %   'Custom', the delay profile is configured using the following 
        %   properties: PathDelays, AveragePathGains, AnglesAoD, AnglesAoA, 
        %   AnglesZoD, AnglesZoA, HasLOSCluster, KFactorFirstCluster, 
        %   AngleSpreads, XPR, NumStrongestClusters.
        %
        %   The default value of this property is 'CDL-A'.
        DelayProfile = 'CDL-A';
        
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
        
        %AnglesAoD Azimuth of departure angles vector (deg)
        %   Specify the azimuth of departure angle for each cluster in
        %   degrees as a double-precision, real, scalar or row vector. This
        %   property applies when DelayProfile is set to 'Custom'.
        %
        %   The default value of this property is 0.0. 
        AnglesAoD = 0.0;
        
        %AnglesAoA Azimuth of arrival angles vector (deg)
        %   Specify the azimuth of arrival angle for each cluster in
        %   degrees as a double-precision, real, scalar or row vector. This
        %   property applies when DelayProfile is set to 'Custom'.
        %
        %   The default value of this property is 0.0.
        AnglesAoA = 0.0;
        
        %AnglesZoD Zenith of departure angles vector (deg)
        %   Specify the zenith of departure angle for each cluster in
        %   degrees as a double-precision, real, scalar or row vector. This
        %   property applies when DelayProfile is set to 'Custom'.
        %
        %   The default value of this property is 0.0.
        AnglesZoD = 0.0;
        
        %AnglesZoA Zenith of arrival angles vector (deg)
        %   Specify the zenith of arrival angle for each cluster in degrees
        %   as a double-precision, real, scalar or row vector. This
        %   property applies when DelayProfile is set to 'Custom'.
        %
        %   The default value of this property is 0.0. 
        AnglesZoA = 0.0;
        
    end
    
    properties (Access = public, Nontunable, Logical)
        
        %HasLOSCluster Line of sight cluster (logical)
        %   Set this property to true to specify that the delay profile has
        %   a line of sight (LOS) cluster. This property applies when
        %   DelayProfile is set to 'Custom'.
        %
        %   The default value of this property is false. 
        HasLOSCluster = false;
        
    end
        
    properties (Access = public, Nontunable)
        
        %KFactorFirstCluster K-factor of first cluster (dB)
        %   Specify the K-factor of the first cluster of the delay profile
        %   in dB (K_1) as a scalar. This property applies when
        %   DelayProfile is set to 'Custom' and HasLOSCluster is set to
        %   true.
        %
        %   The default value of this property is 13.3dB. This is the value
        %   defined for delay profile CDL-D.
        KFactorFirstCluster = 13.3;
        
    end
    
    properties (Access = public, Nontunable, Logical)
        
        %AngleScaling Enable angle scaling (logical)
        %   Set this property to true to apply scaling of angles as
        %   described in TR 38.901 Section 7.7.5.1. This property applies
        %   when you set the DelayProfile property to 'CDL-A', 'CDL-B',
        %   'CDL-C', 'CDL-D' or 'CDL-E'.
        %
        %   The default value of this property is false.
        AngleScaling = false;
        
    end
    
    properties (Access = public, Nontunable)
        
        %AngleSpreads Desired scaled angle spreads vector (deg)
        %   Specify the desired angle spreads (AS_desired) to use for angle
        %   scaling as a row vector [ASD ASA ZSD ZSA]. See TR 38.901
        %   Section 7.7.5.1. ASD is the desired azimuth of departure spread
        %   after scaling. ASA, ZSD and ZSA are the corresponding values 
        %   for azimuth of arrival, zenith of departure and zenith of 
        %   arrival respectively. The values also specify 
        %   [C_ASD C_ASA C_ZSD C_ZSA] for scaling ray offset angles in 
        %   Section 7.7.1 step 1. C_ASD is the desired azimuth angle spread
        %   of departure of rays within a cluster in degrees. C_ASA, C_ZSD
        %   and C_ZSA are the corresponding values for azimuth of arrival, 
        %   zenith of departure and zenith of arrival respectively. This 
        %   property applies when AngleScaling is set to true. Also, when
        %   you set the DelayProfile property to 'Custom' this property is
        %   used to provide the values [C_ASD C_ASA C_ZSD C_ZSA] for
        %   scaling ray offset angles in Section 7.7.1 step 1, but angle
        %   scaling according to Section 7.7.5.1 is not performed.
        %
        %   The default value of this property is [5.0 11.0 3.0 3.0]. This
        %   is the set of values defined for delay profile CDL-A.
        AngleSpreads = [5.0 11.0 3.0 3.0];
        
        %MeanAngles Desired mean scaled angles vector (deg)
        %   Specify the desired mean angles (mu_desired) to use for angle
        %   scaling as a row vector [AoD AoA ZoD ZoA]. See TR 38.901
        %   Section 7.7.5.1. AoD is the desired mean azimuth angle of
        %   departure after scaling. AoA, ZoD and ZoA are the corresponding
        %   values for azimuth of arrival, zenith of departure and zenith
        %   of arrival respectively. This property applies when
        %   AngleScaling is set to true.
        %
        %   The default value of this property is zero for each angle.
        MeanAngles = [0.0 0.0 0.0 0.0];
        
        %XPR Cross polarization power ratio (dB)
        %   Specify the cross-polarization power ratio in dB as a scalar.
        %   This property applies when DelayProfile is set to 'Custom'.
        %
        %   The default value of this property is 10.0dB. This is the value
        %   defined for delay profile CDL-A.
        XPR = 10.0;
        
        %DelaySpread Desired delay spread (s)
        %   Specify the desired RMS delay spread in seconds (DS_desired) as
        %   a scalar. See TR 38.901 Section 7.7.3, and Tables 7.7.3-1 and
        %   7.7.3-2 for examples of desired RMS delay spreads. This
        %   property applies when you set the DelayProfile property to
        %   'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D' or 'CDL-E'.
        %
        %   The default value of this property is 30e-9.
        DelaySpread = 30e-9;
        
        %CarrierFrequency Carrier frequency (Hz)
        %   Specify the carrier frequency in Hertz as a scalar.
        %   
        %   The default value of this property is 4GHz.
        CarrierFrequency = 4e9;
        
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
        
        %UTDirectionOfTravel User terminal (UT) direction of travel
        %   Specify the UT (User Terminal) direction of travel as a column
        %   vector with azimuth and elevation components 
        %   [azimuth; elevation] in degrees.
        %
        %   The default value of this property is [0; 90] degrees.
        UTDirectionOfTravel = [0.0; 90.0];
        
    end
    
    properties (Access = public, Nontunable, Logical)
        
        %KFactorScaling Apply K-factor scaling (logical)
        %   Set this property to true to apply K-factor scaling as
        %   described in TR 38.901 Section 7.7.6. Note that K-factor
        %   scaling modifies both the path delays and path powers. This
        %   property applies if DelayProfile is set to 'CDL-D' or 'CDL-E'.
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
        %   the power of the LOS part of the first cluster to the total
        %   power of all the Laplacian clusters, including the Laplacian
        %   part of the first cluster.
        %
        %   The default value of this property is 9.0dB.
        KFactor = 9.0;
        
        %SampleRate Input signal sample rate (Hz)
        %   Specify the sample rate of the input signal in Hz as a double
        %   precision, real, positive scalar.
        %
        %   The default value of this property is 30.72e6Hz.
        SampleRate = 30.72e6;  
        
        %TransmitAntennaArray Transmit antenna array characteristics
        %   Structure specifying the transmit antenna array. It contains
        %   the following fields:
        %   Size                - Size of antenna array [M,N,P,Mg,Ng]. M
        %                         and N are the number of rows and columns
        %                         in the antenna array. P is the number of
        %                         polarizations (1 or 2). Mg and Ng are the
        %                         number of row and column array panels
        %                         respectively. The defaults are 
        %                         [2,2,2,1,1].
        %   ElementSpacing      - Element spacing in wavelengths expressed
        %                         as [lambda_v lambda_h dg_v dg_h]
        %                         representing the vertical and horizontal
        %                         element spacing and the vertical and
        %                         horizontal panel spacing respectively.
        %                         The defaults are [0.5 0.5 1.0 1.0].
        %   PolarizationAngles  - Polarization angles [theta rho] in
        %                         degrees applicable when P = 2. The
        %                         defaults are [45 -45] degrees.
        %   Orientation         - Array orientation [alpha; beta; gamma] in
        %                         degrees (bearing, downtilt, slant). The
        %                         default values are [0; 0; 0] degrees.
        %   Element             - Antenna element radiation pattern. One of
        %                         'isotropic' or '38.900' (see Section 7.3
        %                         in TR 38.900 and note that the same
        %                         pattern is defined in TR 38.901). The
        %                         default value is '38.900'.
        %   PolarizationModel   - Model describing how to determine the
        %                         radiation field patterns based on a 
        %                         defined radiation power pattern (see
        %                         Section 7.3.2 of TR 38.901). One of
        %                         'Model-1' or 'Model-2'. The default value
        %                         is 'Model-2'.
        %
        % The antenna array elements are mapped to the input waveform
        % channels (columns) in the order that a 5-D array of size
        % M-by-N-by-P-by-Mg-by-Ng is linearly indexed, i.e. across the
        % dimensions first to last. The size of the array is given by
        % TransmitAntennaArray.Size = [M,N,P,Mg,Ng]. For example, an
        % antenna array of size [4,8,2,2,2] will have the first M = 4
        % channels mapped to the first column of the first polarization
        % angle of the first panel. The next M = 4 antennas are mapped to
        % the next column and so on, such that the first M*N = 32 channels
        % are mapped to the first polarization angle of the complete first
        % panel. Then the next 32 channels are mapped in the same fashion
        % to the second polarization angle for the first panel. Subsequent
        % sets of M*N*P = 64 channels are then mapped to the remaining
        % panels, panel rows first then panel columns.
        TransmitAntennaArray = nr5gCDLChannel.antennaArrayStructure([2 2 2 1 1],[0.5 0.5 1.0 1.0],[45 -45],[0;0;0],'38.900',nr5gCDLChannel.defaultPolarizationModel());
        
        %ReceiveAntennaArray Receive antenna array characteristics
        %   Structure specifying the receive antenna array. It contains
        %   the following fields:
        %   Size                - Size of antenna array [M,N,P,Mg,Ng]. M
        %                         and N are the number of rows and columns
        %                         in the antenna array. P is the number of
        %                         polarizations (1 or 2). Mg and Ng are the
        %                         number of row and column array panels
        %                         respectively. The defaults are 
        %                         [1,1,2,1,1].
        %   ElementSpacing      - Element spacing in wavelengths expressed
        %                         as [lambda_v lambda_h dg_v dg_h]
        %                         representing the vertical and horizontal
        %                         element spacing and the vertical and
        %                         horizontal panel spacing respectively.
        %                         The defaults are [0.5 0.5 0.5 0.5].
        %   PolarizationAngles  - Polarization angles [theta rho] in
        %                         degrees applicable when P = 2. The
        %                         defaults are [0 90] degrees.
        %   Orientation         - Array orientation [alpha; beta; gamma] in
        %                         degrees (bearing, downtilt, slant). The
        %                         default values are [0; 0; 0] degrees.
        %   Element             - Antenna element radiation pattern. One of
        %                         'isotropic' or '38.900' (see Section 7.3
        %                         in TR 38.900 and note that the same
        %                         pattern is defined in TR 38.901). The
        %                         default value is 'isotropic'.
        %   PolarizationModel   - Model describing how to determine the
        %                         radiation field patterns based on a 
        %                         defined radiation power pattern (see
        %                         Section 7.3.2 of TR 38.901). One of
        %                         'Model-1' or 'Model-2'. The default value
        %                         is 'Model-2'.
        %
        % The antenna array elements are mapped to the output waveform
        % channels (columns) in the order that a 5-D array of size
        % M-by-N-by-P-by-Mg-by-Ng is linearly indexed, i.e. across the
        % dimensions first to last. The size of the array is given by
        % ReceiveAntennaArray.Size = [M,N,P,Mg,Ng]. For example, an antenna
        % array of size [4,8,2,2,2] will have the first M = 4 channels
        % mapped to the first column of the first polarization angle of the
        % first panel. The next M = 4 antennas are mapped to the next
        % column and so on, such that the first M*N = 32 channels are
        % mapped to the first polarization angle of the complete first
        % panel. Then the next 32 channels are mapped in the same fashion
        % to the second polarization angle for the first panel. Subsequent
        % sets of M*N*P = 64 channels are then mapped to the remaining
        % panels, panel rows first then panel columns.
        ReceiveAntennaArray = nr5gCDLChannel.antennaArrayStructure([1 1 2 1 1],[0.5 0.5 0.5 0.5],[0 90],[0;0;0],'isotropic',nr5gCDLChannel.defaultPolarizationModel());
        
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
        
    end
    
    properties (Access = public, Nontunable)
        
        %SampleDensity Number of time samples per half wavelength
        %   Number of samples of filter coefficient generation per half 
        %   wavelength, i.e. the coefficient generation sampling rate is
        %   F_cg = MaximumDopplerShift * 2 * SampleDensity
        %   Setting SampleDensity = Inf will set 
        %   F_cg = SamplingRate.
        %
        %   The default value of this property is 64.
        SampleDensity = 64;
        
    end
        
    properties (Access = public)
            
        %InitialTime Start time of fading process (s)
        %   Specify the time offset of the fading process as a real
        %   nonnegative scalar. This property is tunable. 
        %
        %   The default value of this property is 0.0.
        InitialTime = 0.0;
        
    end
    
    properties (Access = public, Nontunable)  
        
        %NumStrongestClusters Number of strongest clusters to split into subclusters
        %   The number of strongest clusters to split into sub-clusters.
        %   See TR 38.901 Section 7.5 step 11. This property applies when
        %   DelayProfile is set to 'Custom'.
        %
        %   The default value of this property is 0.
        NumStrongestClusters = 0;
        
        %ClusterDelaySpread Cluster delay spread (s)
        %   Specify the cluster delay spread (C_DS) as a real nonnegative
        %   scalar in seconds. The value is used to specify the delay
        %   offset between sub-clusters for clusters split into
        %   sub-clusters. See TR 38.901 Section 7.5 step 11. This property
        %   applies when DelayProfile is set to 'Custom' and
        %   NumStrongestClusters is greater than zero.
        %
        %   The default value of this property is 3.90625ns.
        ClusterDelaySpread = 5.0/1.28 * 1e-9;
        
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
        
        %ChannelFiltering Perform filtering of input signal (logical)
        %   Set this property to false to disable channel filtering. If set
        %   to false then the step method will not accept an input signal
        %   and the duration of the fading process realization will be
        %   controlled by the NumTimeSamples property (at the sampling rate
        %   given by the SampleRate property). The step method output will
        %   not include an output signal, only the path gains and sample
        %   times.
        %
        %   The default value of this property is true.
        ChannelFiltering = true;
        
    end
        
    properties (Access = public)
        %NumTimeSamples Number of time samples
        %   Specify the number of time samples used to set the duration of
        %   the fading process realization as a positive integer scalar.
        %   This property applies when ChannelFiltering is false. This
        %   property is tunable.
        %
        %   The default value of this property is 30720.
        NumTimeSamples = 30720;
        
    end
        
    properties (Access = public, Nontunable, Logical)
        
        %NormalizeChannelOutputs Normalize channel outputs by the number of receive antennas (logical)
        %   Set this property to true to normalize the channel outputs by
        %   the number of receive antennas. When you set this property to
        %   false, there is no normalization for channel outputs. This
        %   property applies when ChannelFiltering is true.
        %
        %    The default value of this property is true.
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
        
        function set.AnglesAoD(obj,val)
            propName = 'AnglesAoD';
            validateattributes(val,{'double'},{'real','row','finite'},[class(obj) '.' propName], propName);
            obj.AnglesAoD = val;
        end
        
        function set.AnglesAoA(obj,val)
            propName = 'AnglesAoA';
            validateattributes(val,{'double'},{'real','row','finite'},[class(obj) '.' propName], propName);
            obj.AnglesAoA = val;
        end
        
        function set.AnglesZoD(obj,val)
            propName = 'AnglesZoD';
            validateattributes(val,{'double'},{'real','row','finite'},[class(obj) '.' propName], propName);
            obj.AnglesZoD = val;
        end
        
        function set.AnglesZoA(obj,val)
            propName = 'AnglesZoA';
            validateattributes(val,{'double'},{'real','row','finite'},[class(obj) '.' propName], propName);
            obj.AnglesZoA = val;
        end        
        
        function set.KFactorFirstCluster(obj,val)
            propName = 'KFactorFirstCluster';
            validateattributes(val,{'double'},{'real','scalar','finite'},[class(obj) '.' propName],propName);
            obj.KFactorFirstCluster = val;
        end
        
        function set.AngleSpreads(obj,val)
            propName = 'AngleSpreads';
            validateattributes(val,{'double'},{'real','size',[1 4],'finite'},[class(obj) '.' propName], propName);
            obj.AngleSpreads = val;
        end
        
        function set.MeanAngles(obj,val)
            propName = 'MeanAngles';
            validateattributes(val,{'double'},{'real','size',[1 4],'finite'},[class(obj) '.' propName], propName);
            obj.MeanAngles = val;
        end
        
        function set.XPR(obj,val)
            propName = 'XPR';
            validateattributes(val,{'double'},{'real','scalar'},[class(obj) '.' propName],propName);
            obj.XPR = val;
        end
        
        function set.DelaySpread(obj,val)
            propName = 'DelaySpread';
            validateattributes(val,{'double'},{'real','scalar','nonnegative','finite'},[class(obj) '.' propName],propName);
            obj.DelaySpread = val;
        end
        
        function set.CarrierFrequency(obj,val)
            propName = 'CarrierFrequency';
            validateattributes(val,{'double'},{'real','scalar','nonnegative','finite'},[class(obj) '.' propName],propName);
            obj.CarrierFrequency = val;
        end
                
        function set.MaximumDopplerShift(obj,val)
            propName = 'MaximumDopplerShift';
            validateattributes(val,{'double'},{'real','scalar','nonnegative','finite'},[class(obj) '.' propName],propName);
            obj.MaximumDopplerShift = val;
        end
        
        function set.UTDirectionOfTravel(obj,val)
            propName = 'UTDirectionOfTravel';
            validateattributes(val,{'double'},{'real','size',[2 1],'finite'},[class(obj) '.' propName], propName);
            obj.UTDirectionOfTravel = val;
        end
        
        function set.KFactor(obj,val)
            propName = 'KFactor';
            validateattributes(val,{'double'},{'real','scalar','finite'},[class(obj) '.' propName],propName);
            obj.KFactor = val;
        end
        
        function set.TransmitAntennaArray(obj,val)
            obj.TransmitAntennaArray = nr5gCDLChannel.validateAntennaArray(val,class(obj),'TransmitAntennaArray');
        end
        
        function set.ReceiveAntennaArray(obj,val)
            obj.ReceiveAntennaArray = nr5gCDLChannel.validateAntennaArray(val,class(obj),'ReceiveAntennaArray');
        end
        
        function set.SampleRate(obj,val)
            propName = 'SampleRate';
            validateattributes(val,{'double'},{'real','scalar','positive','finite'},[class(obj) '.' propName],propName);
            obj.SampleRate = val;
        end
        
        function set.SampleDensity(obj,val)
            propName = 'SampleDensity';
            validateattributes(val,{'double'},{'real','scalar','positive'},[class(obj) '.' propName],propName);
            obj.SampleDensity = val;
        end
                
        function set.InitialTime(obj,val)
            propName = 'InitialTime';
            validateattributes(val,{'double'},{'real','scalar','nonnegative','finite'},[class(obj) '.' propName],propName);
            obj.InitialTime = val;
        end
        
        function set.NumStrongestClusters(obj,val)
            propName = 'NumStrongestClusters';
            validateattributes(val,{'numeric'},{'scalar','integer','nonnegative'},[class(obj) '.' propName],propName);
            obj.NumStrongestClusters = val;
        end
        
        function set.ClusterDelaySpread(obj,val)
            propName = 'ClusterDelaySpread';
            validateattributes(val,{'double'},{'real','scalar','nonnegative','finite'},[class(obj) '.' propName],propName);
            obj.ClusterDelaySpread = val;
        end
                
        function set.Seed(obj,val)
            propName = 'Seed';
            validateattributes(val,{'double'},{'real','scalar','integer','nonnegative','finite'},[class(obj) '.' propName],propName);
            obj.Seed = val;
        end
        
        function set.NumTimeSamples(obj,val)
            propName = 'NumTimeSamples';
            validateattributes(val,{'numeric'},{'scalar','integer','nonnegative'},[class(obj) '.' propName],propName);
            obj.NumTimeSamples = val;
        end
        
    end
        
    % property value sets for enumerated properties
    properties(Hidden,Transient)
        
        DelayProfileSet = matlab.system.StringSet({'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E','Custom'});
        RandomStreamSet = matlab.system.StringSet({'Global stream','mt19937ar with seed'});
        
    end
    
% =========================================================================
%   protected interface
    
    methods (Access = protected)
        
        % nr5gCDLChannel setupImpl method
        function setupImpl(obj)
            
            % Create underlying channel properties from nr5gCDLChannel
            % properties
            obj.theStruct = nr5gCDLChannel.makeCDLChannelStructure(obj);
            
            % perform parameter validation, perform sub-clustering if
            % applicable, and create information structure
            [obj.theInfo,obj.pdp] = lte.internal.CDLChannelInfo(obj.theStruct);
            
            % Design path delay filters
            if (obj.ChannelFiltering)
                obj.pathFilters = lte.internal.designCDLPathDelayFilters(obj.theInfo.PathDelays,obj.channelFilterDelay,obj.maxFractionalDelayError,obj.SampleRate);
            end
            
            % construct FIRs            
            constructFIRs(obj);

            % setup RNG
            setupRNG(obj);
            
        end
        
        % nr5gCDLChannel stepImpl method
        function varargout = stepImpl(obj,varargin)

            if (obj.ChannelFiltering)
                in = varargin{1};
                insize = size(in);
                outputtype = class(in);
            else
                insize = [obj.NumTimeSamples obj.theInfo.NumTransmitAntennas];
                outputtype = 'double';
            end
            
            [pathgains,sampletimes] = lte.internal.CDLChannel(obj.theStruct,obj.theInfo,obj.pdp,insize,outputtype,obj.rayCoupling,obj.initialPhases,obj.theTime);
            
            if (obj.ChannelFiltering)                
                ntxants = prod(obj.TransmitAntennaArray.Size);
                filtered = applyPathDelayFiltering(obj.FIRs,complex(in(:,1:ntxants)));                                                                   % size 1-by-L cell array of size T-by-P
                out = lte.internal.applyCDLChannelMatrix(obj.SampleDensity,pathgains,sampletimes-obj.theTime,insize,outputtype,obj.SampleRate,filtered); % size T-by-R
                varargout = {out pathgains sampletimes};
            else
                varargout = {pathgains sampletimes};
            end
            
            % advance the time according to the input length, or
            % NumTimeSamples if channel filtering is disabled
            obj.theTime = obj.theTime + (insize(1) / obj.SampleRate);
            
        end
        
        % nr5gCDLChannel resetImpl method
        function resetImpl(obj)
            
            % reset the time to the last InitialTime property value set
            obj.theTime = obj.InitialTime;
            
            % reset filters
            resetFIRs(obj);
            
            % reset RNG
            resetRNG(obj);
            
            % set up initial phases
            obj.initialPhases = lte.internal.getCDLInitialPhases(obj.theStruct,obj.randomStream,obj.theInfo.ClusterTypes);
            
            % set up ray coupling
            obj.rayCoupling = lte.internal.getCDLRayCoupling(obj.theStruct,obj.randomStream,obj.theInfo.ClusterTypes);
            
        end
        
        % nr5gCDLChannel releaseImpl method
        function releaseImpl(obj)
            
            Np = numel(obj.FIRs);
            for np = 1:Np
                release(obj.FIRs{np});
            end            
            obj.pathFilters = [];
            obj.pdp = [];
            obj.initialPhases = [];
            obj.rayCoupling = [];
        
        end
        
        % nr5gCDLChannel getNumInputsImpl method
        function num = getNumInputsImpl(obj)
            
            num = double(obj.ChannelFiltering);
            
        end
        
        % nr5gCDLChannel getNumOutputsImpl method
        function num = getNumOutputsImpl(obj)
            
            num = 2 + obj.ChannelFiltering;
            
        end
        
        % nr5gCDLChannel infoImpl method
        function s = infoImpl(obj)
        %info Returns characteristic information about the CDL channel
        %   S = info(OBJ) returns a structure containing characteristic
        %   information, S, about the CDL fading channel. A description of
        %   the fields and their values is as follows:
        % 
        %   ClusterTypes        - A row cell array of character vectors,
        %                         indicating the type of each cluster in
        %                         the delay profile ('LOS',
        %                         'SubclusteredNLOS', 'NLOS')
        %   PathDelays          - A row vector providing the delays of the
        %                         discrete channel paths, in seconds. These
        %                         values include the effect of the desired
        %                         delay spread scaling, and desired
        %                         K-factor scaling if enabled. 
        %   AveragePathGains    - A row vector of the average gains of the
        %                         discrete path or cluster, in dB. These
        %                         values include the effect of K-factor
        %                         scaling if enabled.
        %   AnglesAoD           - A row vector of the Azimuth of Departure
        %                         angles of the clusters in degrees. These
        %                         values include the effect of angle
        %                         scaling if enabled.
        %   AnglesAoA           - A row vector of the Azimuth of Arrival
        %                         angles of the clusters in degrees. These
        %                         values include the effect of angle
        %                         scaling if enabled.
        %   AnglesZoD           - A row vector of the Zenith of Departure
        %                         angles of the clusters in degrees. These
        %                         values include the effect of angle
        %                         scaling if enabled.
        %   AnglesZoA           - A row vector of the Zenith of Arrival
        %                         angles of the clusters in degrees. These
        %                         values include the effect of angle
        %                         scaling if enabled.
        %   KFactorFirstCluster - K-factor of first cluster of delay
        %                         profile, in dB. If the first cluster of
        %                         the delay profile follows a Laplacian
        %                         rather than Rician distribution,
        %                         KFactorFirstCluster will be -Inf.
        %   NumTransmitAntennas - Number of transmit antennas.
        %   NumReceiveAntennas  - Number of receive antennas.
        %   ChannelFilterDelay  - Channel filter delay in samples.
        %
        %   Note that the step of splitting of the strongest clusters into
        %   sub-clusters described in TR 38.901 Section 7.5 requires
        %   sorting of the clusters by their average power. Therefore if
        %   the NumStrongestClusters property is non-zero (only applicable
        %   for DelayProfile='Custom') the fields of the information
        %   structure are sorted by average power i.e. AveragePathGains is
        %   in descending order of average gain and ClusterTypes,
        %   PathDelays, AnglesAoD, AnglesAoA, AnglesZoD and AnglesZoA are
        %   sorted accordingly. Also, if the HasLOSCluster property is set,
        %   the NLOS (Laplacian) part of that cluster can may be sorted
        %   such that it is not adjacent to the LOS cluster. However,
        %   KFactorFirstCluster will still indicate the appropriate
        %   K-factor.
            
            if isempty(coder.target) && isLocked(obj)
                s = obj.theInfo;
            else
                validate(obj);
                s = lte.internal.CDLChannelInfo(nr5gCDLChannel.makeCDLChannelStructure(obj));
            end
            
        end
        
         % nr5gCDLChannel saveObjectImpl method
        function s = saveObjectImpl(obj)
            
            s = saveObjectImpl@matlab.System(obj);
            s.theStruct = obj.theStruct;
            s.theInfo = obj.theInfo;
            s.pathFilters = obj.pathFilters;
            s.FIRs = obj.FIRs;
            s.randomStream = obj.randomStream;
            s.pdp = obj.pdp;
            s.initialPhases = obj.initialPhases;
            s.rayCoupling = obj.rayCoupling;
            s.theTime = obj.theTime;
            
        end

        % nr5gCDLChannel loadObjectImpl method
        function loadObjectImpl(obj,s,wasLocked)
        
            obj.theTime = s.theTime;
            obj.rayCoupling = s.rayCoupling;
            obj.initialPhases = s.initialPhases;
            obj.pdp = s.pdp;
            obj.randomStream = s.randomStream;
            obj.FIRs = s.FIRs;
            obj.pathFilters = s.pathFilters;
            obj.theInfo = s.theInfo;
            obj.theStruct = s.theStruct;
            loadObjectImpl@matlab.System(obj,s,wasLocked);
            
        end
       
        % nr5gCDLChannel isInactivePropertyImpl method
        function flag = isInactivePropertyImpl(obj,prop)
            
            if (any(strcmp(prop,{'PathDelays','AveragePathGains','AnglesAoD','AnglesAoA','AnglesZoD','AnglesZoA','HasLOSCluster','XPR','NumStrongestClusters'})))
                flag = ~strcmp(obj.DelayProfile,'Custom');
            elseif (strcmp(prop,'KFactorFirstCluster'))
                flag = ~strcmp(obj.DelayProfile,'Custom') || ~obj.HasLOSCluster;
            elseif (strcmp(prop,'AngleSpreads'))
                flag = ~strcmp(obj.DelayProfile,'Custom') && ~obj.AngleScaling;
            elseif (strcmp(prop,'AngleScaling'))
                flag = strcmp(obj.DelayProfile,'Custom');
            elseif (strcmp(prop,'MeanAngles'))
                flag = strcmp(obj.DelayProfile,'Custom') || ~obj.AngleScaling;
            elseif (strcmp(prop,'KFactorScaling'))
                flag = ~any(strcmp(obj.DelayProfile,{'CDL-D','CDL-E'}));
            elseif (strcmp(prop,'KFactor'))
                flag = strcmp(obj.DelayProfile,'Custom') || ~(obj.KFactorScaling);
            elseif (strcmp(prop,'DelaySpread'))
                flag = strcmp(obj.DelayProfile,'Custom');                
            elseif (strcmp(prop,'Seed'))
                flag = ~strcmp(obj.RandomStream,'mt19937ar with seed');
            elseif (strcmp(prop,'ClusterDelaySpread'))
                flag = ~strcmp(obj.DelayProfile,'Custom') || ~obj.NumStrongestClusters;
            elseif (strcmp(prop,'NormalizeChannelOutputs'))
                flag = ~obj.ChannelFiltering;
            elseif (strcmp(prop,'NumTimeSamples'))
                flag = obj.ChannelFiltering;
            else
                flag = false;
            end
            
        end
        
        % nr5gCDLChannel validateInputsImpl method
        function validateInputsImpl(obj,varargin)
            
            if (obj.ChannelFiltering)
                in = varargin{1};
                ntxants = prod(obj.TransmitAntennaArray.Size);
                coder.internal.errorIf(size(in,2)~=ntxants,'lte5g:nr5gCDLChannel:SignalInputNotMatchTxArray',size(in,2),ntxants);
                validateattributes(in,{'double','single'},{'2d','finite'},class(obj),'signal input'); 
            end
            
        end
        
        % nr5gCDLChannel validatePropertiesImpl method
        function validatePropertiesImpl(obj)       
        
            validate(obj);
            
        end

        % nr5gCDLChannel processTunedPropertiesImpl method
        function processTunedPropertiesImpl(obj)
            
            if (isChangedProperty(obj,'InitialTime'))
                
                % if the tuned InitialTime changes the current time (to the
                % same tolerance as used to design fractional delay
                % filters)
                if (abs(obj.theTime - obj.InitialTime) > (obj.maxFractionalDelayError/obj.SampleRate))
                
                    % set the time to the InitialTime property value
                    obj.theTime = obj.InitialTime;

                    % reset filters
                    resetFIRs(obj);
                    
                end
                
            end
            
        end
        
    end
    
    methods (Static, Access = protected)
       
        % validation of antenna array
        function array = validateAntennaArray(val,className,arrayName)
            
            fields = {'Size','ElementSpacing','PolarizationAngles','Orientation','Element'};
            for i = 1:length(fields)
                coder.internal.errorIf(~isfield(val,fields{i}),'lte5g:nr5gCDLChannel:MissingAntennaArrayField',fields{i},arrayName);
            end
            
            classAndArrayName = [className '.' arrayName];
            
            array = struct();
            
            fieldName = 'Size';
            validateattributes(val.Size,{'numeric'},{'integer','size',[1 5],'positive','finite'},[classAndArrayName '.' fieldName],fieldName);
            array.Size = val.Size;
            
            fieldName = 'ElementSpacing';
            validateattributes(val.ElementSpacing,{'double'},{'real','size',[1 4],'nonnegative','finite'},[classAndArrayName '.' fieldName],fieldName);
            array.ElementSpacing = val.ElementSpacing;
            
            fieldName = 'PolarizationAngles';
            validateattributes(val.PolarizationAngles,{'double'},{'real','row','finite'},[classAndArrayName '.' fieldName],fieldName);
            array.PolarizationAngles = val.PolarizationAngles;
            
            fieldName = 'Orientation';
            validateattributes(val.Orientation,{'double'},{'real','size',[3 1],'finite'},[classAndArrayName '.' fieldName],fieldName);
            array.Orientation = val.Orientation;
            
            fieldName = 'Element';
            validatestring(val.Element,{'38.900','isotropic'},[classAndArrayName '.' fieldName],[arrayName '.' fieldName]);
            array.Element = val.Element;
            
            fieldName = 'PolarizationModel';
            if (isfield(val,fieldName)) % New property from R2017b
                validatestring(val.PolarizationModel,{'Model-1','Model-2'},[classAndArrayName '.' fieldName],[arrayName '.' fieldName]);
                array.PolarizationModel = val.PolarizationModel;
            else
                array.PolarizationModel = nr5gCDLChannel.defaultPolarizationModel();
            end
            
        end
        
    end
    
% =========================================================================
%   private

    properties (Access = private, Nontunable)
        
        % The underlying parameter structure used to perform the channel
        % modeling
        theStruct;
        
    end
    
    properties (Access = private)
                
        % Cached info output used to provide info method output
        theInfo;
                
        % path filter coefficients used in channel modeling
        pathFilters;
        
        % path filters
        FIRs;
        
        % RNG stream
        randomStream;
        
        % power delay profile (including columns for angles)
        pdp;
        
        % initial phase array (Phi, TR 38.901 step 10)
        initialPhases;
        
        % ray coupling (TR 38.901 step 8)
        rayCoupling;
        
        % current time, advanced according to the input length (or
        % NumTimeSamples if ChannelFiltering=false) on each step call
        theTime = 0.0;
        
    end
    
    properties (Access = private, Constant)
        
        % channel filter delay, used to design path filters
        channelFilterDelay = 7;
        
        % maximum allowable fractional delay error, used to design path
        % filters and also determine if an InitialTime update requires a
        % filter reset
        maxFractionalDelayError = 0.01;
        
    end
    
    methods (Access = private)
        
        function resetFIRs(obj)
            
            Np = numel(obj.FIRs);
            for np = 1:Np
                reset(obj.FIRs{np});
            end
            
        end
        
        function constructFIRs(obj)
        
            % Note: coder.extrinsic / coder.const is used here in order to
            % make 'coeffs' constant for construction of dsp.FIRFilter
            coder.extrinsic('nr5g.internal.nr5gCDLChannel.makeCDLChannelStructure');
            coder.extrinsic('nr5g.internal.nr5gCDLChannel.getPathDelays');
            coder.extrinsic('lte.internal.designCDLPathDelayFilters');
            if (obj.ChannelFiltering)
                theStruct = coder.const(nr5g.internal.nr5gCDLChannel.makeCDLChannelStructure(obj.NormalizePathGains,obj.NormalizeChannelOutputs,obj.MaximumDopplerShift,obj.UTDirectionOfTravel,obj.CarrierFrequency,obj.Seed,obj.DelaySpread,obj.SampleDensity,obj.SampleRate,obj.DelayProfile,obj.PathDelays,obj.AveragePathGains,obj.AnglesAoD,obj.AnglesAoA,obj.AnglesZoD,obj.AnglesZoA,obj.HasLOSCluster,obj.KFactorFirstCluster,obj.KFactorScaling,obj.KFactor,obj.AngleScaling,obj.AngleSpreads,obj.XPR,obj.ClusterDelaySpread,obj.NumStrongestClusters,obj.MeanAngles,obj.channelFilterDelay,obj.TransmitAntennaArray,obj.ReceiveAntennaArray)); %#ok<PROP>
                pathDelays = coder.const(nr5g.internal.nr5gCDLChannel.getPathDelays(theStruct)); %#ok<PROP>
                [~,coeffs] = coder.const(@lte.internal.designCDLPathDelayFilters,pathDelays,obj.channelFilterDelay,obj.maxFractionalDelayError,obj.SampleRate);
                Np = numel(coeffs);
                obj.FIRs = cell(1,Np);
                for np = 1:Np
                    obj.FIRs{np} = dsp.FIRFilter('Numerator',coeffs{np}.');
                end
            else
                obj.FIRs = {};
            end

        end
        
        function setupRNG(obj)
            
            if strcmp(obj.RandomStream,'Global stream')
                obj.randomStream = [];
            else
                if isempty(coder.target)
                    obj.randomStream = RandStream('mt19937ar','Seed',obj.Seed);
                else
                    obj.randomStream = coder.internal.RandStream('mt19937ar','Seed',obj.Seed);
                end
            end
            
        end
        
        function resetRNG(obj)
            
            if (~strcmp(obj.RandomStream,'Global stream'))
                reset(obj.randomStream,obj.Seed);
            end
            
        end
        
        function validate(obj)
        
            if (strcmp(obj.DelayProfile,'Custom'))
                pdp_props = {'PathDelays','AveragePathGains','AnglesAoD','AnglesAoA','AnglesZoD','AnglesZoA'};
                pdp_rows = validatePDPProperties(obj,pdp_props);
                min_rows = obj.NumStrongestClusters;
                coder.internal.errorIf(pdp_rows < min_rows,'lte5g:nr5gCDLChannel:LessThanStrongestClusters',min_rows,strjoin(pdp_props,','));
            end 
        
        end
        
        function pdp_rows_out = validatePDPProperties(obj,pdp_props)
            
            Nprop = numel(pdp_props);
            pdp_rows = zeros(1,Nprop);
            for i = 1:Nprop
                pdp_rows(i) = numel(obj.(pdp_props{i}));
            end
            coder.internal.errorIf(numel(unique(pdp_rows))~=1,'lte5g:nr5gCDLChannel:UnequalDelayProfileLengths',strjoin(pdp_props,','));
            pdp_rows = unique(pdp_rows);
            pdp_rows_out = pdp_rows(1);
            
        end
        
    end
    
    methods (Static, Access = protected)
        
        % Create underlying channel properties from nr5gCDLChannel
        % properties
        function model = makeCDLChannelStructure(obj)
            
            model = nr5g.internal.nr5gCDLChannel.makeCDLChannelStructure(obj.NormalizePathGains,obj.NormalizeChannelOutputs,obj.MaximumDopplerShift,obj.UTDirectionOfTravel,obj.CarrierFrequency,obj.Seed,obj.DelaySpread,obj.SampleDensity,obj.SampleRate,obj.DelayProfile,obj.PathDelays,obj.AveragePathGains,obj.AnglesAoD,obj.AnglesAoA,obj.AnglesZoD,obj.AnglesZoA,obj.HasLOSCluster,obj.KFactorFirstCluster,obj.KFactorScaling,obj.KFactor,obj.AngleScaling,obj.AngleSpreads,obj.XPR,obj.ClusterDelaySpread,obj.NumStrongestClusters,obj.MeanAngles,obj.channelFilterDelay,obj.TransmitAntennaArray,obj.ReceiveAntennaArray);
            
        end
        
    end
        
    methods (Static, Access = private)
        
        function s = antennaArrayStructure(arraySize,arrayElementSpacing,polarizationAngles,orientation,element,polarizationModel)

            s.Size = arraySize;
            s.ElementSpacing = arrayElementSpacing;
            s.PolarizationAngles = polarizationAngles;
            s.Orientation = orientation;
            s.Element = element;
            s.PolarizationModel = polarizationModel;

        end
        
        function d = defaultPolarizationModel()
            
            d = 'Model-2';
            
        end
    
    end
    
end

function out = applyPathDelayFiltering(firs,in)

    Np = numel(firs);
    if (isempty(in))
        out = repmat({in},1,Np);
    else
        out = cell(1,Np);
        for np = 1:Np
            out{np} = firs{np}(in);
        end
    end
    
end
