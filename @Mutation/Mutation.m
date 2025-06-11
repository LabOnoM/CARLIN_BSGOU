classdef Mutation
    
    properties (SetAccess = immutable, GetAccess = public)
        type
        loc_start
        loc_end
        seq_new
        seq_old
    end
    
    methods (Access = public)
        
        function obj = Mutation(CARLIN_def, type, loc_start, loc_end, seq_new, seq_old)
            
            assert(loc_start >= 1 && loc_start <= CARLIN_def.width.CARLIN);
            assert(type=='M' || type=='I' || type=='D' || type=='C');
            assert(length(seq_old) == length(seq_new));
            
            if (type=='C')
                assert(length(seq_old)>1);
                assert(seq_old(1) ~= seq_new(1) && seq_old(end) ~= seq_new(end));
                seq_old = degap(seq_old);
                seq_new = degap(seq_new);
                max_len = max(length(seq_old), length(seq_new));                
                obj.seq_old = pad(seq_old, max_len, 'right', '-');
                obj.seq_new = pad(seq_new, max_len, 'right', '-');
            else
                if (type=='I')
                    assert(loc_start == loc_end);
                    assert(CARLIN_def.seq.CARLIN(loc_start) == seq_old(1)   && all(seq_old(2:end)  == '-') || ...
                           CARLIN_def.seq.CARLIN(loc_start) == seq_old(end) && all(seq_old(1:end-1)== '-'));
                end
                if (type=='M')
                    assert(loc_start == loc_end);
                    assert(length(seq_old) == 1);
                    assert(CARLIN_def.seq.CARLIN(loc_start) == seq_old);
                    assert(seq_new ~= '-' && seq_new ~= seq_old);
                end
                if (type=='D')
                    assert(isequal(CARLIN_def.seq.CARLIN(loc_start:loc_end), seq_old));
                    assert(all(seq_new=='-'));                
                end            
                obj.seq_old = seq_old;
                obj.seq_new = seq_new;
            end
            obj.type = type;
            obj.loc_start = loc_start;
            obj.loc_end = loc_end;
        end
        
        function s = annotate(obj, full)            
            if (nargin < 2)
                full = false;
            end
            if (obj.type == 'M')
                s = sprintf('%d%s>%s', obj.loc_start, obj.seq_old, obj.seq_new);
            elseif (obj.type == 'D')
                s = sprintf('%d_%ddel', obj.loc_start, obj.loc_end);
            elseif (obj.type == 'I')
                if (full)
                    if (startsWith(obj.seq_old(1), '-'))
                        s = sprintf('%d_%dins%s', obj.loc_start, obj.loc_start+1, degap(obj.seq_new(1:end-1)));
                    else
                        s = sprintf('%d_%dins%s', obj.loc_start, obj.loc_start+1, degap(obj.seq_new(2:end)));
                    end
                else
                    s = sprintf('%d_%dins%d', obj.loc_start, obj.loc_start+1, length(degap(obj.seq_new(2:end))));
                end
            elseif (obj.type == 'C')
                if (full)
                    s = sprintf('%d_%ddelins%s', obj.loc_start, obj.loc_end, degap(obj.seq_new));
                else                    
                    s = sprintf('%d_%ddelins%d', obj.loc_start, obj.loc_end, length(degap(obj.seq_new)));
                end            
            end
        end
        
        function [num_bp_del, num_bp_ins] = num_bps_indel(obj)

            num_bp_del = 0;
            num_bp_ins = 0;

            if (obj.type=='D' || obj.type=='C' || obj.type=='M')
                num_bp_del = (obj.loc_end - obj.loc_start + 1);                
            end
            if (obj.type=='I' || obj.type=='C' || obj.type=='M')        
                num_bp_ins = length(degap(obj.seq_new));
            end
        end
        
    end
    
    methods (Static)
        
        function obj = FromAnnotation(CARLIN_def, str)
           assert(ischar(str));           
           L = length(str);
           
           [obj, s, e] = regexp(str, '(?<loc_start>\d+)(?<seq_old>[ACGT])>(?<seq_new>[ACGT])', 'names');
           if (~isempty(s) && ~isempty(e) && s == 1 && e == L)                              
               obj.loc_start = str2num(obj.loc_start);
               assert(obj.loc_start >= 1 && obj.loc_start <= CARLIN_def.width.CARLIN);
               obj.loc_end = obj.loc_start;
               obj.type = 'M';
               obj = Mutation(CARLIN_def, obj.type, obj.loc_start, obj.loc_end, obj.seq_new, obj.seq_old);
               return;
           end
           
           [obj, s, e] = regexp(str, '(?<loc_start>\d+)_(?<loc_end>\d+)delins(?<seq_new>[ACGT]+)', 'names');           
           if (~isempty(s) && ~isempty(e) && s == 1 && e == L)
               obj.loc_start = str2num(obj.loc_start);
               obj.loc_end = str2num(obj.loc_end);
               assert(obj.loc_start >= 1 && obj.loc_start <= obj.loc_end && obj.loc_end <= CARLIN_def.width.CARLIN);
               obj.seq_old = CARLIN_def.seq.CARLIN(obj.loc_start:obj.loc_end);               
               assert(obj.seq_old(1) ~= obj.seq_new(1) && obj.seq_old(end) ~= obj.seq_new(end));                
               max_len = max(length(obj.seq_old), length(obj.seq_new));                
               obj.seq_old = pad(obj.seq_old, max_len, 'right', '-');
               obj.seq_new = pad(obj.seq_new, max_len, 'right', '-');
               obj.type = 'C';
               obj = Mutation(CARLIN_def, obj.type, obj.loc_start, obj.loc_end, obj.seq_new, obj.seq_old);
               return;
           end
           
           [obj, s, e] = regexp(str, '(?<loc_start>\d+)_(?<loc_end>\d+)ins(?<seq_new>[ACGT]+)', 'names');
           if (~isempty(s) && ~isempty(e) && s == 1 && e == L)
               obj.loc_start = str2num(obj.loc_start);
               obj.loc_end = str2num(obj.loc_end);
               assert(obj.loc_start >= 0 && obj.loc_end==obj.loc_start+1 && obj.loc_end <= CARLIN_def.width.CARLIN+1);
               insert_left = obj.loc_start==0 || ismembc(obj.loc_start, CARLIN_def.bounds.consites(:,2));
               if (insert_left)
                   assert(obj.loc_end <= CARLIN_def.width.CARLIN);
                   obj.loc_start = obj.loc_end;
                   obj.seq_new = [obj.seq_new CARLIN_def.seq.CARLIN(obj.loc_end)];
                   obj.seq_old = pad(CARLIN_def.seq.CARLIN(obj.loc_end), length(obj.seq_new), 'left', '-');
               else
                   assert(obj.loc_start >= 1);
                   obj.loc_end = obj.loc_start;
                   obj.seq_new = [CARLIN_def.seq.CARLIN(obj.loc_start) obj.seq_new];
                   obj.seq_old = pad(CARLIN_def.seq.CARLIN(obj.loc_start), length(obj.seq_new), 'right', '-');
               end
               obj.type = 'I';
               obj = Mutation(CARLIN_def, obj.type, obj.loc_start, obj.loc_end, obj.seq_new, obj.seq_old);
               return;
           end
           
           [obj, s, e] = regexp(str, '(?<loc_start>\d+)_(?<loc_end>\d+)del', 'names');
           if (~isempty(s) && ~isempty(e) && s == 1 && e == L)
               obj.loc_start = str2num(obj.loc_start);
               obj.loc_end = str2num(obj.loc_end);
               assert(obj.loc_start >= 1 && obj.loc_start <= obj.loc_end && obj.loc_end <= CARLIN_def.width.CARLIN);
               obj.seq_old = CARLIN_def.seq.CARLIN(obj.loc_start:obj.loc_end);               
               obj.seq_new = repmat('-', [1 length(obj.seq_old)]);
               obj.type = 'D';
               obj = Mutation(CARLIN_def, obj.type, obj.loc_start, obj.loc_end, obj.seq_new, obj.seq_old);
               return;
           end
           error('Invalid mutation annotation: %s', str);
            
        end
        
        function mut_list = FromFile(CARLIN_def, filename)
            assert(exist(filename, 'file') == 2, sprintf('Missing annotation file: %s', filename));
            s = splitlines(fileread(filename));
            if (isempty(s{end}))
                s = s(1:end-1);
            end
            mut_list = cell(size(s));
            for i = 1:length(mut_list)                
                if (strcmp(s{i},'?') || strcmp(s{i},'X'))
                    mut_list{i} = s{i};
                elseif (~strcmp(s{i}, '[]'))                    
                    mut_list{i} = cellfun(@(x) Mutation.FromAnnotation(CARLIN_def, x), strsplit(s{i}, ',')', 'un', false);
                    mut_list{i} = vertcat(mut_list{i}{:});
                end
            end
        end

        function ToFile(CARLIN_def, mut_list, folder, filename)
            
            if (~exist(folder, 'dir'))
                mkdir(folder);        
            end
            
            fid = fopen([folder '/' filename], 'wt');
            
            if (~isempty(mut_list))
                if (isa(mut_list{1}, 'AlignedSEQ'))
                    mut_list = cellfun(@(x) Mutation.identify_cas9_events(CARLIN_def, x), mut_list, 'un', false);                
                end
                s = cellfun(@(x) arrayfun(@(i) x(i).annotate(true), [1:size(x,1)], 'un', false), mut_list, 'un', false);
                s = cellfun(@(x) strjoin(x,','), s, 'un', false);
                s(cellfun(@isempty, s)) = {'[]'};

                fprintf(fid, '%s\n', s{:});
            end
            
            fclose(fid);
            
        end
        
        m = make_compound_mutation(CARLIN_def, s, e, in_seq, ref_seq);
        modified_sites = find_modified_sites(CARLIN_def, in);
        [mut_list, bp_event, seq, ref_seq, ref_mask] = identify_sequence_events(CARLIN_def, aligned);
        [event_list, mutated_allele] = identify_cas9_events(CARLIN_def, aligned);
        mutated_allele = apply_mutations(CARLIN_def, event_list);
        [bp_event, del_event, ins_event] = classify_bp_event(CARLIN_def, aligned);
        
        [num_bps_del_mu, num_bps_ins_mu, num_bps_del_se, num_bps_ins_se] = num_bps_indel_stats(summary, edited_only);
        
    end
end