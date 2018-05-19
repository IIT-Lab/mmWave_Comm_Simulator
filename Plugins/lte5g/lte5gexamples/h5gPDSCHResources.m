%h5gPDSCHResources PDSCH and DM-RS resource element indices and DM-RS values
%   [IND,DMRSIND,DMRS,INFO] = h5gPDSCHResources(GNB,CHS) returns the 
%   resource element (RE) indices for a PDSCH transmission and associated
%   PDSCH DM-RS, given the time (symbols) and frequency (PRBs) allocation
%   of the PDSCH and the DM-RS configuration. The 1-based linear PDSCH
%   indices are returned in matrix, IND. They are defined relative to a 
%   three-dimensional RE grid representing a 14-symbol slot for the full
%   carrier (in the PDSCH numerology) across the layers/DM-RS ports of the
%   PDSCH. Each column of IND represents the grid locations for a separate 
%   layer/port (the third dimension of the grid). The DM-RS RE indices 
%   have the same format and are returned in matrix, DMRSIND. The complex 
%   values of DM-RS sequence are also returned in matrix, DMRS. 
%  
%   The cell-wide setting input, GNB, must be a structure including
%   the fields:
%   NDLRB  - Number of downlink resource blocks for the carrier 
%            (in the PDSCH numerology)
%
%   The PDSCH specific input, CHS, must be a structure including the
%   fields:
%   NSlot               - Slot number of PDSCH transmission
%   PRBSet              - PRBs allocated to the PDSCH (0-based indices)
%   SymbolSet           - OFDM symbols allocated to the PDSCH within the slot
%                         (including DM-RS symbols, 0-based indices, max range 0...13)
%   PortSet             - DM-RS ports used by PDSCH
%                         (0-based indices, max range 0...11,
%                         mapping to ports p=1000...1011 respectively)
%   NLayers             - Number of layers
%   Modulation          - Modulation ('QPSK','16QAM','64QAM','256QAM')
%   Reserved            - Reserved PRB patterns (structure array, see below)
%   DL_DMRS_config_type - DM-RS configuration type (1,2)
%   NIDNSCID            - DM-RS scrambling identity (0...65535)
%   NSCID               - DM-RS scrambling initialization (0,1)
% 
%   The DM-RS symbol locations can either be defined explicitly using the
%   single parameter:
%   DMRSSymbolSet       - OFDM symbols containing the DM-RS within the PDSCH
%                         allocation (0-based indices, max range 0...13)
%   or, defined implicitly via the following group of parameters:
%   PDSCHMappingType    - PDSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
%   DL_DMRS_typeA_pos   - Mapping type A only. First DM-RS symbol position (2,3)
%   DL_DMRS_max_len     - Number of front-loaded DM-RS symbols 
%                         (1(single symbol),2(double symbol))
%   DL_DMRS_add_pos     - Additional DM-RS symbol positions (max range 0...3)
%
%   Periodically recurring patterns of reserved PRB can defined using the 
%   'Reserved' parameter. These PRB will be excluded from the generated 
%   indices and the DL-SCH/PDSCH processing should rate-match around them. 
%   It can be used to exclude SS/PBCH, CORESETs/PDCCH and other resources,
%   as defined in TS 38.214 section 5.1.4. This parameter takes the format
%   of a structure array where each element defines a separate pattern. 
%   Each element, i, of the array should contain the following fields:
%   Reserved(i).PRB     - Reserved PRB (0-based indices, defined as a 
%                         vector or cell array)
%   Reserved(i).Symbols - OFDM symbols associated with reserved PRB 
%                         (0-based indices, spanning one or more slots)
%   Reserved(i).Period  - Total number of slots in the pattern period
% 
%   The reserved PRB indices can be specified as a vector or a cell array.
%   If a vector then the same PRBs are excluded in each OFDM symbol in the 
%   pattern. If the PRB indices are defined as a cell array then each cell
%   specifies the excluded PRBs for the associated OFDM symbol in the 
%   pattern. In the latter case, the length of the PRB cell array should
%   match the length of the 'Symbols' field, i.e. an individual set of PRBs
%   is defined for each reserved symbol. The symbols that form the 
%   time-domain locations of the reserved pattern can be greater than 13
%   and therefore cover multiple slots. The overall pattern will repeat
%   itself every 'Period' slots. If this field is empty then the pattern 
%   will cyclically repeat itself after the whole number of slots covered 
%   by the specified symbols. 
% 
%   In terms of frequency domain DM-RS density, there are two different RRC
%   signaled configuration types ('DL_DMRS_config_type'). Configuration
%   type 1 defines 6 subcarriers per PRB per antenna port, comprising
%   alternating subcarriers. Configuration type 2 defines 4 subcarriers per
%   PRB per antenna ports, consisting of 2 groups of 2 neighbouring
%   subcarriers. Different shifts are applied to the sets of subcarriers
%   used, dependent on the associated antenna port or CDM group. For
%   type 1, there are 2 possible CDM groups/shifts across up to 8 possible
%   antenna ports (p=1000...1007), and, for type 2, there are 3 possible
%   CDM groups/shifts across 12 ports (p=1000...1011). See TS 38.211
%   section 7.4.1.1 for the full configuration details.
%
%   In terms of the time-domain DM-RS symbol positions, the PDSCH mapping
%   type ('PDSCHMappingType') can be either slot-wise (type A) or non
%   slot-wise (type B). When a UE is scheduled to receive PDSCH by a DCI,
%   this mapping type is signaled by the time-domain resource field in the
%   grant. The field acts as an index into an RRC configured table where
%   each row in the table specifies a combination of mapping type, slot
%   offset, K0, the symbol start and length indicator, SLIV. The mapping
%   type specifies the relative locations of the associated DM-RS. For
%   slot-wise mapping type A, the first DM-RS symbol is signaled by a field
%   in the MIB to be either 2 or 3 ('DL_DMRS_typeA_pos'). For the non
%   slot-wise mapping type B, the first DM-RS symbol is always the first
%   symbol of the PDSCH time allocation.
% 
%   The maximum number of DM-RS OFDM symbols used by a UE is configured by
%   RRC signalling ('DL_DMRS_add_pos' and 'DMRS_max_len'). The DM-RS can be 
%   a set of single symbols, distributed roughly uniformly across the 
%   allocated PDSCH symbols, or 1 or 2 pairs of neighbouring or 'double 
%   symbol' DM-RS. The 'DMRS_max_len' RRC parameter (1 or 2 respectively)
%   configures whether only single symbol DM-RS or either single or double
%   symbol DM-RS are used. In the latter case, the actual selection is 
%   signaled in the DCI format 1_1 message. The 'DL_DMRS_add_pos'
%   higher-layer parameter defines the number of single or double symbol
%   DM-RS that are transmitted. The valid combinations of these two 
%   parameters is given by TS 38.211 tables 6.4.1.1.3-3 and 6.4.1.1.3-4.
%   In this function, the value of the 'DL_DMRS_max_len' input parameter
%   directly controls whether either single or double symbols are used. 
% 
%   INFO is the output structure containing the fields:
%   G             - Bit capacity of the PDSCH. This is should be the
%                   length of codeword to be output from the DL-SCH 
%                   transport channel
%   Gd            - Number of resource elements per layer/port, equal to 
%                   the number of rows in the PDSCH indices
%   DMRSSymbolSet - The symbol numbers in a slot containing DM-RS (0-based)
%   NREPerPRB     - Number of RE per PRB allocated to PDSCH (not accounting
%                   for any reserved resources)
%
%   Example:
%   % Display the locations of the PDSCH and PDSCH DM-RS resource elements 
%   % in a slot. 
%     
%   % Set the number of downlink carrier resource blocks (in the 
%   % PDSCH numerology) 
%   gnb = struct('NDLRB',50); 
% 
%   % Specific the basic PDSCH allocation properties to be full band, full
%   % slot using 2 layers/ports
%   pdsch = struct();
%   pdsch.NSlot = 0;               % Slot number
%   pdsch.PRBSet = 0:gnb.NDLRB-1;  % Full band PRBs allocated to the PDSCH
%   pdsch.SymbolSet = 0:13;        % Full slot allocation
%   pdsch.PortSet = [0,2];         % Use DM-RS ports 1000 and 1002
%   pdsch.Modulation = 'QPSK';     % Modulation type
%
%   % Exclude some PRBs (e.g. for a CORESET) in the middle of the carrier
%   % in the first two symbols in the slot
%   pdsch.Reserved.PRB = fix(gnb.NDLRB/2) + (-4:5);
%   pdsch.Reserved.Symbols = [0,1];
%   pdsch.Reserved.Period = 1;
% 
%   % Configure the PDSCH DM-RS for config type 1 and slot-wise, type A
%   % mapping. Use double symbols for front-loaded DM-RS and configure an
%   % additional symbol pair towards the end of the slot
%   pdsch.DL_DMRS_config_type = 1; % DM-RS configuration type (1,2)
%   pdsch.NIDNSCID = 1;            % DM-RS scrambling identity (0...65535)
%   pdsch.NSCID = 0;               % DM-RS scrambling initialization (0,1)
%   pdsch.PDSCHMappingType = 'A';  % Slot-wise PDSCH mapping type
%   pdsch.DL_DMRS_typeA_pos = 2;   % First DM-RS symbol position for type A 
%   pdsch.DL_DMRS_max_len = 2;     % Specify double front-loaded DM-RS
%   pdsch.DL_DMRS_add_pos = 1;     % Specify an additional DM-RS pair
% 
%   % Display PDSCH and DM-RS RE locations of the first port of the grid
%   slotgrid = zeros(12*gnb.NDLRB,14,length(pdsch.PortSet));
%   [ind,dmrsind,dmrs,info] = h5gPDSCHResources(gnb,pdsch);
%   slotgrid(ind) = 20;                 % Use light blue for PDSCH RE 
%   slotgrid(dmrsind) = 50*abs(dmrs);   % Use yellow for DM-RS RE
%   figure;
%   image(slotgrid(:,:,1));
%   title('PDSCH and DM-RS resource elements');
%   axis('xy'), xlabel('Symbols'), ylabel('Subcarriers');
%  
%   See also h5gDLSCH, ltePDSCHIndices, lteDMRSIndices.

%   Copyright 2018 The MathWorks, Inc.

function [pdschIndices,dmrsIndices,dmrssymb,pdschInfo] = h5gPDSCHResources(gnb, pdsch)
   
    % Check the input argument number
    narginchk(2,2);

    % Capture the set of OFDM symbols (within the slot) associated with PDSCH 
    if isfield(pdsch,'SymbolSet')
        allocatedsymbols = pdsch.SymbolSet;
    else
        allocatedsymbols = (0:13);              % Full slot if not specified
    end
           
    % Capture the set of OFDM symbols (within the slot) that will carry the DM-RS
    % If not explicitly defined then calculate DM-RS containing OFDM symbols from parameter set
    if isfield(pdsch,'DMRSSymbolSet')
        dmrssymbols = intersect(pdsch.DMRSSymbolSet,allocatedsymbols);
    else
        dmrssymbols = getDMRSymbols(pdsch,allocatedsymbols);
    end
    
    % Capture the set of antenna ports associated with the PDSCH and DM-RS
    % Config type 1: p=1000-1007, Config type 2: p=1000-1011 (max ranges)
    if isfield(pdsch,'PortSet')
        ports = reshape(pdsch.PortSet,1,[]);
    else
        ports = 0:pdsch.NLayers-1;
    end    
    
    % Capture the slot number of the PDSCH
    % Used to locate the current slot in relation to any extended reserved resource patterns
    if isfield(pdsch,'NSlot')
       nslot = pdsch.NSlot;
    else
       nslot = 0;
    end
    
    % RE level mapping (PDSCH containing RE per PRB)
    % Find all PDSCH RE per PRB for each OFDM symbol that will also carry DM-RS
    % 
    configtype1 = pdsch.DL_DMRS_config_type==1;
    if (configtype1)
        dmrssc = [0 2 4 6 8 10]';        % Configuration type 1: 6 DM-RS symbols per PRB per CDM 
        dshifts = mod(fix(ports/2),2);   % Delta shift per CDM group (0,1)... the *set* of shifts/groups that we will exclude
    else      
        dmrssc = [0 1 6 7]';             % Configuration type 2: 4 DM-RS symbols per PRB per CDM
        dshifts = 2*mod(fix(ports/2),3); % Delta shift per CDM group (0,2,4)... the *set* of shifts/groups that we will exclude    
    end
    fullprb = ones(12,1);        % Map all the subcarriers in a PRB
    dmrsre = dmrssc + dshifts;   % Implicit expansion of zero-shifted RE positions across all the shifts/CDM in play
    fullprb(dmrsre+1) = 0;       % Clear all RE which will carry DM-RS in at least one port
    pdschre = find(fullprb)-1;   % Find PDSCH (non DM-RS) RE in a DM-RS containing symbol
       
    % PRB level mapping (PDSCH containing PRB)
    % Find all allocated PDSCH PRB (accounting for any reserved PRB) across
    % the allocated OFDM symbols
    %  
    prbset = reshape(pdsch.PRBSet,1,[]);   % Ensure that it is a row for implicit expansion later
    prbcell = cell(1,14);
    prbcell(allocatedsymbols+1) = {prbset};
    
    % Loop over the reserved resources and exclude them from the PDSCH PRB
    % across the allocated symbols
    if isfield(pdsch,'Reserved')
        for ri=1:length(pdsch.Reserved)

            % Reference the current reserved symbol/PRB indices pair and period
            reservedsymbols = pdsch.Reserved(ri).Symbols;
            reservedprb = reshape(pdsch.Reserved(ri).PRB,1,[]);
            reservedperiod = pdsch.Reserved(ri).Period;

            % Find any of the allocated symbols which overlap with reservations
            %
            % If the reserved period was empty then get number of complete slots
            % in reservation period and cyclically extend pattern to cover current slot 
            if isempty(reservedperiod)
                reservedperiod = ceil((max(reservedsymbols)+1)/14);      % Period, in terms of whole slots
            end
            offset = mod(nslot,reservedperiod)*14; 
            inter = intersect(allocatedsymbols,reservedsymbols-offset);  % Find allocated symbols which contain reserved PRB
            % Reference the PRB associated with the overlapping symbols
            if iscell(reservedprb)
                if length(reservedprb) ~= length(reservedsymbols)  
                    error('Symbol-wise reserved PRB resources (defined in a cell array) and associated symbols should have the same number of elements.');
                end
                prbdiff = arrayfun(@(x)setdiff(prbcell{x+1},[reservedprb{x==(reservedsymbols-offset)}]), inter,'UniformOutput',false);
            else
                prbdiff = arrayfun(@(x)setdiff(prbcell{x+1},reservedprb), inter,'UniformOutput',false);
            end
            prbcell(inter+1) = prbdiff;
        end   
    end
        
    % Create the RE per PRB list across the allocated symbols, accounting
    % for DM-RS and non DM-RS carrying cases
    recell = cell(1,14);
    recell(allocatedsymbols+1) = {(0:11)'};  % PDSCH RE per PRB in normal data symbols
    recell(dmrssymbols+1) = {pdschre};       % PDSCH RE per PRB PRB in DM-RS carrying symbols
        
    % Combine PRB oriented and RE per PRB oriented index arrays and expand (using
    % implicit expansion) into a column of linear indices for a single antenna/layer
    % Implicit expansion:                         Column             Row                                      Column
    slotindices = cell2mat( arrayfun(@(x) reshape(recell{x+1} + 12*prbcell{x+1} + 12*gnb.NDLRB*x,[],1), allocatedsymbols(:),'UniformOutput',false) );
     
    % Expand across all antenna ports/layers and make 1-based
    pdschIndices = slotindices(:) + (1 + (12*14*gnb.NDLRB)*(0:numel(ports)-1));
    
    % PDSCH bit capacity information
    pdschInfo.Gd = size(pdschIndices,1);         % Number of QAM symbols in one layer/antenna grid
    q = 2*find(strcmpi(pdsch.Modulation,{'QPSK','16QAM','64QAM','256QAM'}));
    pdschInfo.G = q*pdschInfo.Gd*length(ports);  % And scale by layers/ports and modulation scheme
    pdschInfo.DMRSSymbolSet = dmrssymbols;       % Symbols containing DM-RS
    pdschInfo.NREPerPRB  = 12*(length(allocatedsymbols)-length(dmrssymbols)) + length(dmrssymbols)*length(pdschre);

    % Create the DM-RS QPSK symbols and resource element indices associated
    % with the PDSCH transmission
    %
    % First create unshifted version on the DM-RS RE indices in a single antenna plane
    % Implicit expansion:                           Column        Row                                  Column
    dmslotindices = cell2mat( arrayfun(@(x) reshape(dmrssc + 12*prbcell{x+1} + 12*gnb.NDLRB*x,[],1), dmrssymbols(:),'UniformOutput',false) );
    
    % Expand across all antenna ports/layers, applying DM-RS shifts for each port, and make 1-based
    % Implicit expansion: Column             Row                             Row
    dmrsIndices = dmslotindices(:) + (1 + dshifts + (12*14*gnb.NDLRB)*(0:numel(ports)-1));
    
    % Create the matching set of DM-RS QPSK symbols for all ports
    dmrssymb = getDMRSSequence(pdsch,prbset,dmrssymbols);
     
end

% Get DM-RS complex symbol values
function symb = getDMRSSequence(pdsch,prbset,dmrssymbols)
    
    nidnscid = pdsch.NIDNSCID;  % DL-DMRS-Scrambling-ID
    nscid = pdsch.NSCID;        
    
    if isfield(pdsch,'PortSet')
        ports = pdsch.PortSet;
    else 
        ports = 0:pdsch.NLayers-1;
    end
       
    % Capture the slot number of the PDSCH
    % Used to initialize the scrambling cinit and to locate the current slot
    % in relation to any extended reserved resource patterns
    if isfield(pdsch,'NSlot')
       nslot = pdsch.NSlot;
    else
       nslot = 0;
    end

    % Deduce whether double DM-RS symbols are in use
    doublesymbol =  length(dmrssymbols)>1 && (dmrssymbols(end)-dmrssymbols(end-1))==1 && mod(length(dmrssymbols),2)==0;

    % DM-RS config type (1 or 2)
    type1 = pdsch.DL_DMRS_config_type==1;
    nsc = 8 + type1*4;      % Number of PRBS bits for 6 DM-RS QPSK symbols (type 1) or 4 DM-RS QPSK symbols (type 2) per PRB
    prbcap = max(prbset)+1; % Find the maximum PRB indices for the allocation

    % We will cache the 'outer' set of PRBS in this array 
    prbs = zeros(length(prbset)*nsc,length(dmrssymbols));

    % Expand the PDSCH PRB indices into 1-based 
    % Implicit expansion: Column     Row
    scidx = reshape((1:nsc)' + nsc*prbset(:)',[],1);  

    % PRB level mapping (PDSCH containing PRB)
    % Find all allocated PDSCH PRB (accounting for any reserved PRB) across
    % the allocated OFDM symbols
    prbidx = cell(1,size(prbs,2));
    % Loop over DM-RS containing symbols
    for i=1:size(prbs,2)

        % Generate the PRBS sequence associated with PRB = 0...max(PRB)
        % 5G uses the same PRBS generator as LTE
        cinit = mod(2^17*(14*mod(nslot,10) + dmrssymbols(i) + 1)*(2*nidnscid + 1) + 2*nidnscid + nscid,2^31);
        prbscol = ltePRBS(cinit,nsc*prbcap);    
        % Extract only the PRBS portions associated with the allocated PRB set
        prbs(:,i) = prbscol(scidx);

        % Remove any reserved PRB which overlaps with allocation 
        rprb = [];
        % Loop over the reserved resources
        if isfield(pdsch,'Reserved')
            for ri=1:length(pdsch.Reserved)

                % Reference the current reserved symbol/PRB indices pair and period
                reservedsymbols = pdsch.Reserved(ri).Symbols;
                reservedprb = reshape(pdsch.Reserved(ri).PRB,1,[]);
                reservedperiod = pdsch.Reserved(ri).Period;

                % Find any of the allocated symbols which overlap with reservations
                %
                % If the reserved period was empty then get number of complete slots
                % in reservation period and cyclically extend pattern to cover current slot 
                if isempty(reservedperiod)
                    reservedperiod = ceil((max(reservedsymbols)+1)/14);      % Period, in terms of whole slots
                end
                offset = mod(nslot,reservedperiod)*14;   
                % Logical vector corresponding any reserved PRB for current DM-RS symbol
                rsymb = dmrssymbols(i)==(reservedsymbols-offset);  

                % Combine any reserved PRB for this DM-RS symbol into a single array
                if any(rsymb)
                  if iscell(reservedprb)
                    rprb = [rprb reservedprb{rsymb}]; %#ok<AGROW> Get all reserved PRB values associated with current DM-RS symbol
                  else
                    rprb = [rprb reservedprb];        %#ok<AGROW>
                  end
                end  
            end   
        end
        
        % Remove any reserved PRB from the allocated PRB for the current DM-RS symbol
        if ~isempty(rprb)
            [~,prbidx{i}] = setdiff(prbset,rprb,'stable');
        else
            prbidx{i} = (1:length(prbset))';
        end
    end
    
    % Modulate extracted portions (CMD multiplications could also be implemented as logical operations)
    % Reshape into a rectangular block (subcarriers-by-symbols)
    basesymb = reshape(lteSymbolModulate(prbs(:),'QPSK'),size(prbs,1)/2,size(prbs,2));

    % Turn the DM-RS symbol wise indices, relative to PRB allocation (accounting for any reservations),
    % into full 1-based linear indices associated with extracted DM-RS symbol values
    % Implicit expansion:                           Column                Row
    slotindices = cell2mat( arrayfun(@(x) reshape((1:nsc/2)' + (nsc/2)*(prbidx{x+1}-1)' + size(basesymb,1)*x,[],1),(0:size(basesymb,2)-1)','UniformOutput',false) );

    % Apply time/frequency masks
    %
    % Frequency mask is applied to every other DM-RS symbol in frequency
    % A time mask is only applied to every other double symbols
    % 
    % Allocate array for the output
    symb = zeros(numel(slotindices),numel(ports));

    % Loop over the ports and apply masks
    for pidx = 1:numel(ports)
        
        pv = ports(pidx);
        % Mask applied if p is odd or (p > 4 or 6 and double-symbol)
        if mod(pv,2)
          m1 = [1; -1];
        else
          m1 = 1;      
        end
        if doublesymbol && pv > 6
          m2 = [1 -1];
        else
          m2 = 1;
        end
        m = m1*m2;

        % Apply mask to QPSK symbols
        dmrs = basesymb.*repmat(m,size(basesymb)./size(m));

        % Extract the port-dependent DM-RS symbol values associated 
        % with the used PRB for the current slot
        symb(:,pidx) = dmrs(slotindices);
    end

end

% Get the OFDM symbol numbers contain DM-RS for the PDSCH allocation
function dmrssymbols = getDMRSymbols(pdsch,allocatedsymbols)        
    
    % Get PDSCH mapping type (A (slot-wise) or B (non slot-wise))
    typeB = strcmpi(pdsch.PDSCHMappingType,'B');
    
    % Establish the first DM-RS symbol according to the mapping type
    if typeB 
        dmrssymbols = allocatedsymbols(1);     % If type B (non slot-wise)       
    else
        dmrssymbols = pdsch.DL_DMRS_typeA_pos; % If type A (slot-wise) then 2 or 3
    end
     
    % Double-symbol DM-RS (in the operating system, RRC configured with length 2 then DCI signaled) 
    if pdsch.DL_DMRS_max_len == 2
        dmrssymbols = [dmrssymbols, dmrssymbols+1];
    end
   
    % Create static tables for mapping type (A/B) and single or double-symbol   
    persistent dmrs_add_pos;
    if isempty(dmrs_add_pos) 
    
        % Additional position tables
        % Single-symbol, 1,2,3 *additional* symbols
        dmrs_singleA = {
         [],[],    [];
         [],[],    [];
         [],[],    [];
         7, [],    [];
         9, [6,9], [];
         9, [6,9], [];
         9, [6,9], [5,8,11];
         11,[7,11],[5,8,11];
         11,[7,11],[5,8,11];
        };
        dmrs_singleB = {
         [],[],[];
         [],4 ,[];
         [],[],[];
         [],[],[];
         [],[],[];
         [],[],[];
         [],[],[];
         [],[],[];
         [],[],[];
        }; 
        % Double-symbol, 1,2 *additional* symbol pairs
        dmrs_doubleA = {
         [],     [];
         [],     [];
         [],     [];
         [],     [];
         [8,9],  [];
         [8,9],  [];
         [8,9],  [];
         [10,11],[];
         [10,11],[];
        };
        dmrs_doubleB = {
         [],[];
         [],[];
         [],[];
         [],[];
         [],[];
         [],[];
         [],[];
         [],[];
         [],[];
        };   

        % Combined tables, indexed as tables{type,length}
        dmrs_add_pos = { dmrs_singleA, dmrs_doubleA; dmrs_singleB, dmrs_doubleB };  
    end
    
    % Add any additional DM-RS symbol positions as defined
    if pdsch.DL_DMRS_add_pos
        addpos = dmrs_add_pos{ typeB+1, pdsch.DL_DMRS_max_len }; 
        if pdsch.DL_DMRS_add_pos <= size(addpos,2)
            nsymbols = length(allocatedsymbols);
            tix = max(1,nsymbols-5);  % 1-based index wrt PDSCH duration into additional DM-RS positions table
            dmrssymbols = [dmrssymbols, addpos{tix,pdsch.DL_DMRS_add_pos}];
        end
    end
    
end
