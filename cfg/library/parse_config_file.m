function cfg = parse_config_file(cfg_file)

    preconfig = ~endsWith(cfg_file, '.json');
    
    if (preconfig)
        [folder, ~, ~] = fileparts(mfilename('fullpath'));    
        fprintf('Parsing CFG file: %s.json\n', cfg_file);

        assert(ismember(cfg_file, {'Sanger'; 'BulkDNA'; 'BulkRNA'; 'sc10xV2'; 'sc10xV3'; 'scInDropsV2'; 'scInDropsV3'}), ...
               'Unrecognized CFG Type. Choose from: {Sanger, BulkDNA, BulkRNA, sc10xV2, sc10xV3, scInDropsV2, scInDropsV3}');

        cfg_file = sprintf('%s/%s.json', folder, cfg_file);
    else        
        fprintf('Parsing CFG file: %s\n', cfg_file);
    end
    
    assert(exist(cfg_file, 'file') == 2, sprintf('Missing config file: %s', cfg_file));
    cfg = jsondecode(fileread(cfg_file));
    
    assert(isfield(cfg, 'type'), 'Missing type field');
    assert(ismember(cfg.type, {'Sanger'; 'Bulk'; 'SC'}), 'Valid values for type field: {Sanger, Bulk, SC}');
    
    assert(isfield(cfg, 'UMI') && isfield(cfg.UMI, 'length'), 'Missing UMI.length field');
    assert(isscalar(cfg.UMI.length), 'UMI.length should be an integer value');
    
    if (isequal(cfg.type, 'SC'))
        
        assert(isfield(cfg, 'SC'), 'Missing SC field');
        assert(isfield(cfg.SC, 'Platform'), 'Missing SC.Platform field');
        assert(isfield(cfg.SC, 'Version'), 'Missing SC.Version field');
        
        assert(ismember(cfg.SC.Platform, {'InDrops'; '10x'}), 'Valid values for SC.Platform: {InDrops, 10x}');
        assert(ismember(cfg.SC.Version, [2; 3]), 'Valid values for SC.Version: {2, 3}');
        
        assert(isfield(cfg, 'CB'), 'Missing CB field');
        assert(isfield(cfg.CB(1), 'length'), 'Missing CB(1).length field');
        assert(isnumeric(cfg.CB(1).length), 'CB(1).length should be integer-valued');
        
        assert(~isfield(cfg.UMI, 'location'), 'cfg.UMI.location should not be included if cfg.type=SC');
        
        if (isequal(cfg.SC.Platform, 'InDrops'))
            assert(length(cfg.CB)==2, 'For InDrops, CB field should be an array with 2 elements');
            assert(isfield(cfg.CB(2), 'length'), 'Missing CB(2).length field');
            assert(isnumeric(cfg.CB(2).length), 'CB(2).length should be integer-valued');
        end
        
    else
        assert(isfield(cfg.UMI, 'location'), 'Missing UMI.location field');
        assert(ismember(cfg.UMI.location, {'L'; 'R'}), 'Valid values for UMI.location field: {L, R}');
        assert(~isfield(cfg, 'SC'), 'SC field should not be set when cfg.type=%s', cfg.type);
        assert(~isfield(cfg, 'CB'), 'CB field should not be specified cfg.type=%s', cfg.type);
    end
    
    assert(isfield(cfg, 'trim'), 'Missing trim field');
    assert(isfield(cfg.trim, 'Primer5'), 'Missing trim.Primer5 field');
    assert(isfield(cfg.trim, 'Primer3'), 'Missing trim.Primer3 field');
    assert(isfield(cfg.trim, 'SecondarySequence'), 'Missing trim.SecondarySequence field');
    
    assert(ismember(cfg.trim.Primer5, {'exact'; 'misplaced'; 'malformed'; 'ignore'}), ...
        'Valid values for trim.Primer5 field: {exact, misplaced, malformed, ignore}');
    assert(ismember(cfg.trim.Primer3, {'exact'; 'misplaced'; 'malformed'; 'ignore'}), ...
        'Valid values for trim.Primer3 field: {exact, misplaced, malformed, ignore}');
    assert(ismember(cfg.trim.SecondarySequence, {'exact'; 'misplaced'; 'malformed'; 'ignore'}), ...
        'Valid values for trim.SecondarySequence field: {exact, misplaced, malformed, ignore}');
    
    assert(isfield(cfg, 'read_perspective'), 'Missing read_perspective field');
    assert(isfield(cfg.read_perspective, 'ShouldComplement'), 'Missing read_perspective.ShouldComplement field');
    assert(isfield(cfg.read_perspective, 'ShouldReverse'), 'Missing read_perspective.ShouldReverse field');
    
    assert(ismember(cfg.read_perspective.ShouldComplement, {'N'; 'Y'}), ...
        'Valid values for read_perspective.ShouldComplement field: {N, Y}');
    assert(ismember(cfg.read_perspective.ShouldReverse, {'N'; 'Y'}), ...
        'Valid values for read_perspective.ShouldReverse field: {N, Y}');
  
end
