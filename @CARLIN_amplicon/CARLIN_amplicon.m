classdef (Sealed) CARLIN_amplicon < handle
        
    %% Constant Definitions
    properties (SetAccess = immutable)
        N;
        motifs;
        seq;
        width;
        bounds;
        bps;        
        match_score;
        open_penalty;
        close_penalty;
        sub_score;
    end 
    
    methods (Access = public)
        
        function obj = CARLIN_amplicon(con)
            
            assert(nargin == 1, 'CARLIN_amplicon needs to called with a CARLIN amplicon configuration file');
            
            % A lot of redundant definitions here, but helpful to compute
            % them up front once in the singleton definition, rather than
            % repeatedly on the fly later.
            
            % Set counts of different element types
            obj.N.segments = size(con.sequence.segments,1);
            obj.N.pams     = obj.N.segments-1;
            obj.N.motifs   = 2*obj.N.segments + obj.N.pams + 2;
            
            % Define types of different motifs
            obj.motifs.prefix   = 1;
            obj.motifs.consites = [2:3:obj.N.motifs];
            obj.motifs.cutsites = [3:3:obj.N.motifs];
            obj.motifs.pams     = [4:3:obj.N.motifs-1];
            obj.motifs.postfix  = obj.N.motifs;
            
            % Define standard sequences
            obj.seq.Primer5 = con.sequence.Primer5;
            obj.seq.Primer3 = con.sequence.Primer3;
            obj.seq.SecondarySequence = con.sequence.SecondarySequence;            
            
            obj.seq.prefix   = con.sequence.prefix;
            obj.seq.postfix  = con.sequence.postfix;
            obj.seq.pam      = con.sequence.pam;            
            obj.seq.segments = con.sequence.segments;
            temp = cell(obj.N.segments+obj.N.pams+2,1);
            temp(1) = cellstr(obj.seq.prefix);
            temp(2:2:end) = obj.seq.segments;
            temp(3:2:end-1) = obj.seq.pam;
            temp(end) = cellstr(obj.seq.postfix);
            obj.seq.CARLIN = horzcat(temp{:});
            
            % Precompute widths
            obj.width.Primer5 = length(obj.seq.Primer5);
            obj.width.Primer3 = length(obj.seq.Primer3);
            obj.width.SecondarySequence = length(obj.seq.SecondarySequence);
            obj.width.prefix  = length(obj.seq.prefix);
            obj.width.segment = length(obj.seq.segments{1});
            obj.width.postfix = length(obj.seq.postfix);
            if (isempty(obj.seq.pam))
                obj.width.pam = 0;
            else
                obj.width.pam = length(obj.seq.pam{1});
            end
            obj.width.CARLIN  = length(obj.seq.CARLIN);
            obj.width.consite = size(con.open_penalty.consites,2);
            obj.width.cutsite = size(con.open_penalty.cutsites,2);            
            obj.width.min_length = obj.width.prefix+obj.width.consite+obj.width.postfix;
            
            % Define custom sequences
            obj.seq.consites = cellfun(@(x) x(1:obj.width.consite), obj.seq.segments, 'un', false);
            obj.seq.cutsites = cellfun(@(x) x(obj.width.consite+1:end), obj.seq.segments, 'un', false);
           
            % Define endpoints of different elements
            temp = cumsum(cellfun(@length, temp));
            temp = [[1; 1+temp(1:end-1)] temp];
            
            obj.bounds.prefix   = temp(1,:);
            obj.bounds.segments = temp(2:2:end,:);
            obj.bounds.pams     = temp(3:2:end-1,:);
            obj.bounds.postfix  = temp(end,:);
        
            obj.bounds.consites  = [obj.bounds.segments(:,1) obj.bounds.segments(:,2)-obj.width.cutsite];
            obj.bounds.cutsites  = [obj.bounds.segments(:,1)+obj.width.segment-obj.width.cutsite obj.bounds.segments(:,2)];            
            
            % Precompute bounds of ordered motifs
            obj.bounds.ordered                        = zeros(obj.N.motifs,2);
            obj.bounds.ordered(obj.motifs.prefix,:)   = obj.bounds.prefix;
            obj.bounds.ordered(obj.motifs.consites,:) = obj.bounds.consites;
            obj.bounds.ordered(obj.motifs.cutsites,:) = obj.bounds.cutsites;
            obj.bounds.ordered(obj.motifs.pams,:)     = obj.bounds.pams;
            obj.bounds.ordered(obj.motifs.postfix,:)  = obj.bounds.postfix;
            
            % Precompute BP membership of all motif types
            obj.bps.prefix = obj.bounds.prefix(1):obj.bounds.prefix(2);
            obj.bps.consite = arrayfun(@(s,e) s:e, obj.bounds.consites(:,1)', obj.bounds.consites(:,2)', 'un', false);
            obj.bps.consite = horzcat(obj.bps.consite{:});
            obj.bps.cutsite = arrayfun(@(s,e) s:e, obj.bounds.cutsites(:,1)', obj.bounds.cutsites(:,2)', 'un', false);
            obj.bps.cutsite = horzcat(obj.bps.cutsite{:});
            obj.bps.pam     = arrayfun(@(s,e) s:e, obj.bounds.pams(:,1)', obj.bounds.pams(:,2)', 'un', false);
            obj.bps.pam     = horzcat(obj.bps.pam{:});
            obj.bps.postfix = obj.bounds.postfix(1):obj.bounds.postfix(2);
            
            % Empirically derived NUC44 alignment score thresholds to determine
            % a successful match.
            obj.match_score.Primer5   = con.match_score.Primer5;
            obj.match_score.Primer3   = con.match_score.Primer3;
            obj.match_score.SecondarySequence = con.match_score.SecondarySequence;
            
            open_penalty = cell(obj.N.motifs,1);
            open_penalty(obj.motifs.prefix )  = {con.open_penalty.prefix};
            open_penalty(obj.motifs.postfix)  = {con.open_penalty.postfix};
            open_penalty(obj.motifs.cutsites) = num2cell(con.open_penalty.cutsites, 2);
            open_penalty(obj.motifs.consites) = num2cell(con.open_penalty.consites, 2);
            open_penalty(obj.motifs.pams)     = num2cell(con.open_penalty.pam, 2);            
            obj.open_penalty = [con.open_penalty.init horzcat(open_penalty{:})];
    
            close_penalty = cell(obj.N.motifs,1);
            close_penalty(obj.motifs.prefix )  = {con.close_penalty.prefix};
            close_penalty(obj.motifs.postfix)  = {con.close_penalty.postfix};
            close_penalty(obj.motifs.cutsites) = num2cell(con.close_penalty.cutsites, 2);
            close_penalty(obj.motifs.consites) = num2cell(con.close_penalty.consites, 2);
            close_penalty(obj.motifs.pams)     = num2cell(con.close_penalty.pam, 2);
            obj.close_penalty = horzcat(close_penalty{:});

            sub_score = nuc44;
            obj.sub_score = sub_score(1:4,1:4);
            
        end
           
        % Helper functions
        loc = locate(ref, bp, motif_bounds);                
        motif = find_motif_from_bp(obj, n);
        site = coarse_grain_motif(obj, n);

        % Deconstruct aligned sequences
        chop_sites = partition_sequence(obj, seq_aligned, ref_aligned);
        [out, chop_sites] = desemble_sequence(obj, seq_aligned, ref_aligned);        

        % Alignment
        [sc, al] = cas9_align(obj, seq);        
     
    end
        
    methods (Static)
        
        % Helper functions
        ordered_motifs = order_named_motifs(ref, motifs);        
        
        % Alignment
        [sc, al] = cas9_align_mex(seq, ref, open_penalty, close_penalty, sub_score);
        [so1, so2, ro] = pairwise_align(a1, a2);        
        sc = nwalign_score(al, match_score, mismatch_penalty, gap_open, gap_extend);
        
        % Construct aligned sequence
        out = assemble_sequence(segments);
        
    end
end
