function RGB = get_sequence_coloring(CARLIN_def, aligned_seqs, resolution)

    assert(iscell(aligned_seqs) && isa(aligned_seqs{1}, 'AlignedSEQ'));
    if (nargin < 2)
        resolution = 'bp';
    end
    
    switch(resolution)
        case 'bp'
            event = cellfun(@(x) Mutation.classify_bp_event(CARLIN_def, x), aligned_seqs, 'Un', false);
            event = vertcat(event{:});
        case 'motif'
            event = cellfun(@(x) AlignedSEQMotif.classify_motif_event(x.get_seq(), x.get_ref()), aligned_seqs, 'Un', false);
            event = vertcat(event{:});
            event = repelem(event, 1, diff(CARLIN_def.bounds.ordered,[],2)+1);
        otherwise
            error("Unrecognized sequence coloring resolution");
    end
    
    RGB = arrayfun(@(x) CARLIN_viz.color(x), event, 'un', false);
    RGB = reshape(vertcat(RGB{:}), [], CARLIN_def.width.CARLIN, 3);
    
end
