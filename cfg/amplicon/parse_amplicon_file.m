function con = parse_amplicon_file(amplicon_file)

    if(~endsWith(amplicon_file, '.json'))
        fprintf('Parsing CARLIN amplicon file: %s.json\n', amplicon_file);
        [folder, ~, ~] = fileparts(mfilename('fullpath'));
        if (ismember(amplicon_file, {'OriginalCARLIN'; 'TigreCARLIN'}))
            amplicon_file = sprintf('%s/%s.json', folder, amplicon_file);
        else
            error('Unrecognized preset amplicon file: %s. Valid options are {OriginalCARLIN, TigreCARLIN}', amplicon_file)
        end
    else
        fprintf('Parsing CARLIN amplicon file: %s\n', amplicon_file);
    end
    
    assert(exist(amplicon_file, 'file') == 2, 'Missing CARLIN amplicon file: %s', amplicon_file);    
  
    con = jsondecode(fileread(amplicon_file));
    
    fields = {'sequence'; 'match_score'; 'open_penalty'; 'close_penalty'};
    
    for i = 1:length(fields)    
        assert(isfield(con, fields{i}), 'Missing %s field', fields{i});
    end
    
    primer_motifs = {'Primer5'; 'Primer3'; 'SecondarySequence'};
    
    for i = 1:length(primer_motifs)
        assert(isfield(con.sequence, primer_motifs{i}), 'Missing sequence.%s field', primer_motifs{i});
        assert(ischar(con.sequence.(primer_motifs{i})), 'sequence.%s should be a character array', primer_motifs{i});
        assert(isfield(con.match_score, primer_motifs{i}), 'Missing match_score.%s field', primer_motifs{i});
        assert(isscalar(con.match_score.(primer_motifs{i})), 'match_score.%s should be scalar', primer_motifs{i});
    end
    
    outer_motifs = {'prefix'; 'postfix'};   
    
    for i = 1:length(outer_motifs)
        
        assert(isfield(con.sequence, outer_motifs{i}), 'Missing sequence.%s field', outer_motifs{i});
        assert(ischar(con.sequence.(outer_motifs{i})), 'sequence.%s should be a character array', outer_motifs{i});
        
        assert(isfield(con.open_penalty, outer_motifs{i}), 'Missing open_penalty.%s field', outer_motifs{i});
        assert(isnumeric(con.open_penalty.(outer_motifs{i})), 'open_penalty.%s should be numeric', outer_motifs{i});
        assert(length(con.sequence.(outer_motifs{i})) == length(con.open_penalty.(outer_motifs{i})), ...
            'open_penalty.%s needs to be the same length as sequence.%s', outer_motifs{i}, outer_motifs{i});
        con.open_penalty.(outer_motifs{i}) = con.open_penalty.(outer_motifs{i})';
        
        assert(isfield(con.close_penalty, outer_motifs{i}), 'Missing close_penalty.%s field', outer_motifs{i});
        assert(isnumeric(con.close_penalty.(outer_motifs{i})), 'close_penalty.%s should be numeric', outer_motifs{i});        
        assert(length(con.sequence.(outer_motifs{i})) == length(con.close_penalty.(outer_motifs{i})), ...
            'close_penalty.%s needs to be the same length as sequence.%s', outer_motifs{i}, outer_motifs{i});
        con.close_penalty.(outer_motifs{i}) = con.close_penalty.(outer_motifs{i})';
        
    end
       
    assert(isfield(con.sequence, 'segments'), 'Missing sequence.segments field');    
    N_segment = size(con.sequence.segments,1);    
    assert(N_segment>0, 'Must specify at least one target site in sequence.segment');
    if (N_segment==1)         
        if (ischar(con.sequence.segments))
            con.sequence.segments = {con.sequence.segments};
        end
        L_segment = length(con.sequence.segments{1});
    else        
        assert(all(cellfun(@ischar, con.sequence.segments)), 'sequence.segments should be an array of nucleotide strings');
        L_segment = unique(cellfun(@length, con.sequence.segments));        
        assert(isscalar(L_segment), 'All target sites in sequence.segments should have the same length');
    end
    
    assert(isfield(con.sequence, 'pam'), 'Missing sequence.pam field');
    N_pam = size(con.sequence.pam,1);
    if (N_pam==0)
        con.sequence.pam = {};
        assert(N_segment==1, 'Must specify non-empty nucleotide string for sequence.pam when more than one sequence.segment is present');
        L_pam = 0;
    elseif (N_pam==1)
        assert(N_segment>1, 'sequence.pam should be an empty string if sequence.segments specifies a single target site');
        if (ischar(con.sequence.pam))
            con.sequence.pam = {con.sequence.pam};
        end
        con.sequence.pam = repmat(con.sequence.pam, [N_segment-1, 1]);
        L_pam = length(con.sequence.pam{1});
    else
        assert(N_pam==N_segment-1, 'When multiple nucleotide strings are specified for sequence.pam, the # PAM sequences needs to be 1 less than the # of target sites');
        assert(all(cellfun(@ischar, con.sequence.pam)), 'sequence.pam should be an array of nucleotide strings');
        L_pam = unique(cellfun(@length, con.sequence.pam));        
        assert(isscalar(L_pam), 'nucleotide strings in sequence.pam should all be the same length');
    end
    
    [con, L_open_con ] = check_penalty_vector(con, 'open_penalty',  'consites', N_segment);
    [con, L_close_con] = check_penalty_vector(con, 'close_penalty', 'consites', N_segment);    
    assert(L_open_con==L_close_con, 'Length of open_penalty.consites must equal length of close_penalty.consites');
            
    [con, L_open_cut ] = check_penalty_vector(con, 'open_penalty',  'cutsites', N_segment);
    [con, L_close_cut] = check_penalty_vector(con, 'close_penalty', 'cutsites', N_segment);
    assert(L_open_cut==L_close_cut, 'Length of open_penalty.cutsites must equal length of close_penalty.cutsites');
    
    assert(L_open_con+L_open_cut==L_segment, 'Combined length of open_penalty.consites and open_penalty.cutsites must equal length of target site');
    assert(L_close_con+L_close_cut==L_segment, 'Combined length of close_penalty.consites and close_penalty.cutsites must equal length of target site');
    
    [con, L_open_pam ] = check_penalty_vector(con, 'open_penalty',  'pam', N_segment-1);
    [con, L_close_pam] = check_penalty_vector(con, 'close_penalty', 'pam', N_segment-1);
    assert(L_open_pam==L_close_pam, 'Length of open_penalty.pam must equal length of close_penalty.pam');
    
    assert(L_open_pam==L_pam, 'Length of (open_penalty/close_penalty).pam must equal length of sequence.pam');
    
    assert(isfield(con.open_penalty, 'init'), 'Missing open_penalty.init field');
    assert(isscalar(con.open_penalty.init), 'open_penalty.init should be scalar');
    
end

function [con, L_penalty, N_penalty] = check_penalty_vector(con, ptype, site, N_ref)
    assert(isfield(con.(ptype), site), 'Missing %s.%s field', ptype, site);
    assert(isnumeric(con.(ptype).(site)), '%s.%s should be numeric and all rows should have the same length', ptype, site);
    N_penalty = size(con.(ptype).(site));
    if (N_penalty(2) > 1)
        L_penalty = N_penalty(2);
        N_penalty = N_penalty(1);
        if (strcmp(site, 'pam'))
            assert(N_penalty == N_ref, '%s.%s should be a single row or have one less row than the number of target sites in sequence.segments', ptype, site);
        else
            assert(N_penalty == N_ref, '%s.%s should be a single row or have the same number of rows as target sites in sequence.segments', ptype, site);
        end
    else
        L_penalty = N_penalty(1);
        N_penalty = 1;
        con.(ptype).(site) = con.(ptype).(site)';
        con.(ptype).(site) = repmat(con.(ptype).(site), [N_ref, 1]);
    end
    
end
