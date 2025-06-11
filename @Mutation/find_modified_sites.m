function modified_sites = find_modified_sites(CARLIN_def, in)

    assert(isa(in, 'AlignedSEQ'));

    modified_sites = [];
    event_list = Mutation.identify_cas9_events(CARLIN_def, in);
    
    if (~isempty(event_list))
        s_motif = arrayfun(@(i) CARLIN_def.coarse_grain_motif(CARLIN_def.find_motif_from_bp(i)), vertcat(event_list.loc_start));
        e_motif = arrayfun(@(i) CARLIN_def.coarse_grain_motif(CARLIN_def.find_motif_from_bp(i)), vertcat(event_list.loc_end  ));

        modified_sites = arrayfun(@(i) s_motif(i):e_motif(i), 1:length(event_list), 'un', false);
        modified_sites = unique(horzcat(modified_sites{:}));
    end
    
end