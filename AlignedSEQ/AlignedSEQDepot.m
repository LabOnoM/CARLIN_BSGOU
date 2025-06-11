classdef AlignedSEQDepot < handle

    properties (SetAccess = private, GetAccess = public)
        unaligned_SEQ
        aligned_SEQ
        alignment_map
        backing_file = '';
    end
    
    methods (Access = public)
        
        function obj = AlignedSEQDepot(SEQ, CARLIN_def)
            
            assert(~any(cellfun(@(x) any(x=='-'), SEQ)), 'No insertions should be present in sequences to align');
            
            assert(length(SEQ) == length(unique(SEQ)), 'Duplicate sequences should not be supplied to AlignedSEQDepot');            
                        
            N = size(SEQ,1);
            fprintf('Creating aligned CARLIN store with %d sequences\n', N);
            v = cell(N,1);
            parfor i = 1:N
                [~, v{i}] = CARLIN_def.cas9_align(SEQ{i});
            end            
            
            obj.unaligned_SEQ = SEQ;
            obj.aligned_SEQ = v;
            obj.alignment_map = uint32([1:N]');
        end
        
        function alignment = get_alignment_for_SEQ_ind(obj, SEQ_ind)
            if isempty(obj.aligned_SEQ) && ~isempty(obj.backing_file)
                obj.load_from_disk();
            end
            assert(SEQ_ind > 0 && SEQ_ind <= length(obj.alignment_map));
            alignment = obj.aligned_SEQ{obj.alignment_map(SEQ_ind)};
        end
        
        function alignment = get_alignment_for_SEQ(obj, SEQ)
            if isempty(obj.aligned_SEQ) && ~isempty(obj.backing_file)
                obj.load_from_disk();
            end
            is = find(ismember(obj.unaligned_SEQ, SEQ));
            if (isempty(is))
                alignment = [];
            else
                assert(length(is) == 1);
                alignment = obj.aligned_SEQ{obj.alignment_map(is)};
            end
        end
        
        function sanitize_prefix_postfix(obj)
            fprintf('Sanitizing prefix/postfix of aligned sequences in depot\n');            
            v = AlignedSEQ.sanitize_prefix_postfix(obj.aligned_SEQ);
            [~, backtrack, ind] = unique(cellfun(@(x) degap(x.get_seq()), v, 'un', false));
            obj.aligned_SEQ = v(backtrack);
            obj.alignment_map = ind(obj.alignment_map);
            fprintf('...reduced sequence diversity from %d to %d\n', length(ind), length(obj.aligned_SEQ));
        end
        
        function sanitize_conserved_regions(obj, CARLIN_ref)
            fprintf('Sanitizing conserved regions of aligned sequences in depot\n');
            v = AlignedSEQ.sanitize_conserved_regions(obj.aligned_SEQ, CARLIN_ref);
            [~, backtrack, ind] = unique(cellfun(@(x) degap(x.get_seq()), v, 'un', false));
            obj.aligned_SEQ = v(backtrack);
            obj.alignment_map = ind(obj.alignment_map);
            fprintf('...reduced sequence diversity from %d to %d\n', length(ind), length(obj.aligned_SEQ));
        end

        function spill_to_disk(obj, filename)
            if nargin < 2
                filename = [tempname '.mat'];
            end
            s = struct('unaligned_SEQ', {obj.unaligned_SEQ}, 'aligned_SEQ', {obj.aligned_SEQ}, ...
                       'alignment_map', obj.alignment_map);
            save(filename, '-struct', 's', '-v7.3');
            obj.unaligned_SEQ = [];
            obj.aligned_SEQ = [];
            obj.alignment_map = [];
            obj.backing_file = filename;
        end

        function load_from_disk(obj)
            if isempty(obj.backing_file) || exist(obj.backing_file, 'file') ~= 2
                return;
            end
            s = load(obj.backing_file);
            obj.unaligned_SEQ = s.unaligned_SEQ;
            obj.aligned_SEQ = s.aligned_SEQ;
            obj.alignment_map = s.alignment_map;
        end

        function seqs = get_aligned_SEQs(obj)
            if isempty(obj.aligned_SEQ) && ~isempty(obj.backing_file)
                obj.load_from_disk();
            end
            seqs = obj.aligned_SEQ;
        end
        
    end
end