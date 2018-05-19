%h5gCRCEncode Cyclic redundancy check calculation and appending
%   BLKCRC = h5gCRCEncode(...) calculates a CRC for the input data vector
%   and returns a copy of the vector with the CRC attached. As an option
%   the CRC can be masked before appending. To support the processing of
%   filler bits, negative input bit values are interpreted as logical 0 for
%   the purposes of the CRC calculation (-1 is used to represent filler
%   bits).
%  
%   BLKCRC = h5gCRCEncode(BLK,POLY) calculates the CRC defined by POLY for
%   input bit vector BLK and returns a copy of the input with the CRC
%   appended in vector BLKCRC. The CRC polynomial is defined by a value
%   from the set ('6','11','16','24A','24B','24C'). See TS 38.212 Section
%   5.1 for details the associated polynomials.
% 
%   BLKCRC = h5gCRCEncode(BLK,POLY,MASK) behaves as above except the third
%   parameter allows the appended CRC bits to be xor masked with the
%   integral value of MASK. The MASK value is applied to the CRC bits MSB
%   first/LSB last.
%    
%   Example 1: 
%   % The CRC associated with an all zero vector return an all
%   % zero vector of length 124.
%
%   crc1 = h5gCRCEncode(zeros(100,1),'24A');
%   crc1(1:5)
%   
%   Example 2: 
%   % The CRC bits are masked in a MSB first order, resulting
%   % in all zeros apart from a single one in the last element 
%   % position.
%
%   crc2 = h5gCRCEncode(zeros(100,1),'24A',1);
%   crc2(end-5:end)
%   
%   See also h5gCRCDecode, h5gCodeBlockSegment, lteCRCEncode.

%   Copyright 2018 The MathWorks, Inc.

function blkcrc = h5gCRCEncode(blk, poly, mask)

    persistent polylist;
    persistent crclenlist;
    persistent polynomlist;
    persistent encoders;
    
    narginchk(2, 3);
    if nargin == 2
        mask = 0;
    end
    
    % Create encoder on demand
    if isempty(encoders)
        encoders = cell(1,6);
        polynomlist = {[1 1 0 0 0 0 1],...                                     % 6 CRC bits for polar uplink, 18<=K<=25
                       [1 1 1 0 0 0 1 0 0 0 0 1], ...                          % 11 CRC bits for polar uplink, K>30
                       [1 0 0 0 1 zeros(1, 6) 1 0 0 0 0 1], ...                % 16 CRC bits for LDPC transport block size <= 3824
                       [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1],... % 24 CRC bits for LDPC transport block size > 3824
                       [1 1 zeros(1, 16) 1 1 0 0 0 1 1],...                    % 24 CRC bits for LDPC code block segments
                       [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1]};   % 24 CRC bits for polar downlink (BCH and DCI)
        polylist = {'6', '11', '16', '24A', '24B', '24C'};
        crclenlist = {6, 11, 16, 24, 24, 24};
    end

    % Validate input information bits
    if iscell(blk)
        error('lte:error','MATLAB array should be numeric and cannot be of cell type.');
    end
    if isstruct(blk)
        error('lte:error','MATLAB array should be numeric and cannot be of struct type.');
    end
    if ~isvector(blk) && ~isempty(blk)
        error('lte:error','The input must be a vector and not a matrix.');
    end
    blk = blk(:); % comm.CRCGenerator can only handle a column
    
    % Validate CRC type
    encoderInd = strcmpi(poly,polylist);
    if ~any(encoderInd)
        error('lte:error',['The input CRC polynormial type is not one of the set (' strjoin(polylist,', ') ')']);
    end
    
    % Validate mask
    if iscell(mask)
        error('lte:error','MATLAB array should be numeric and cannot be of cell type.');
    end
    if isstruct(mask)
        error('lte:error','MATLAB array should be numeric and cannot be of struct type.');
    end
    if ischar(mask) || isstring(mask)
        error('lte:error','MATLAB array cannot be of char or string type.');
    end
    
    % Set mask to 0 when it is empty
    if isempty(mask)
        mask = 0;
    end
    
    % uint32 conversion is to handle inf or nan inputs,
    % double conversion is for bitshift operation later
    mask = double(uint32(real(mask(1))));
    
    % Create encoder object on demand and perform encoding
    encoder = encoders{encoderInd};
    if isempty(encoder)
       encoder = comm.CRCGenerator('Polynomial',polynomlist{encoderInd});
       encoders{encoderInd} = encoder;
    end
    % Any values larger than 0 are treated as a logical 1, any values
    % smaller or equal to 0 (including filler bits) are treated as a
    % logical 0
    outBits = encoder(blk>0);
    
    % Apply mask to CRC bits and concatenate to the input block
    if isempty(blk)
        blkcrc = zeros(0,1,class(blk));
    else
        crcBitsLen = crclenlist{encoderInd}; 
        maskBits = mod(bitshift(mask,(1-crcBitsLen:0)'),2); % Convert decimal mask to bits
        blkcrc = [blk; xor(outBits(end-crcBitsLen+1:end),maskBits)];
    end

end