function tests = CARLIN_def_test
    tests = functiontests(localfunctions);
end

function test_cutsite_width(testCase)
    check_cutsite_width(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_cutsite_width(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')));
end

function check_cutsite_width(testCase, CARLIN_def)    
    verifyEqual(testCase, CARLIN_def.width.cutsite, 7);
end

function test_CARLIN_length(testCase)
    check_CARLIN_length(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_CARLIN_length(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')));
end

function check_CARLIN_length(testCase, CARLIN_def)    
    verifyEqual(testCase, length(CARLIN_def.seq.CARLIN), 276);
    verifyEqual(testCase, CARLIN_def.width.CARLIN, 276);
end

function test_start_pos(testCase)
    check_start_pos(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_start_pos(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')))
end

function check_start_pos(testCase, CARLIN_def)
    verifyEqual(testCase, CARLIN_def.bounds.prefix(1), 1);
end

function test_end_pos(testCase)
    check_end_pos(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_end_pos(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')));
end

function check_end_pos(testCase, CARLIN_def)    
    verifyEqual(testCase, CARLIN_def.bounds.postfix(2), CARLIN_def.width.CARLIN);
end

function test_continuous_indices_segments(testCase)
    check_continuous_indices_segments(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_continuous_indices_segments(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')));
end

function check_continuous_indices_segments(testCase, CARLIN_def)    
    actual = cell(CARLIN_def.N.segments+CARLIN_def.N.pams+2,1);
    actual{1}         = CARLIN_def.bounds.prefix;
    actual(2:2:end)   = num2cell(CARLIN_def.bounds.segments,2);
    actual(3:2:end-1) = num2cell(CARLIN_def.bounds.pams,2);
    actual{end}       = CARLIN_def.bounds.postfix;
    actual = cellfun(@(x) x(1):x(2), actual, 'un', false);
    actual = horzcat(actual{:});
    verifyEqual(testCase, actual, [1:CARLIN_def.width.CARLIN]);
end

function test_bounds_prefix(testCase)
    check_bounds_prefix(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_bounds_prefix(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')));
end

function check_bounds_prefix(testCase, CARLIN_def)    
    verifyEqual(testCase, CARLIN_def.seq.CARLIN(CARLIN_def.bounds.prefix(1):CARLIN_def.bounds.prefix(2)), CARLIN_def.seq.prefix);
end

function test_bounds_consites(testCase)
    check_bounds_consites(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_bounds_consites(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')));
end

function check_bounds_consites(testCase, CARLIN_def)    
    expected = cellfun(@(x) x(1:CARLIN_def.width.consite), CARLIN_def.seq.segments, 'un', false);    
    actual_from_bounds = arrayfun(@(s,e) CARLIN_def.seq.CARLIN(s:e), CARLIN_def.bounds.consites(:,1), CARLIN_def.bounds.consites(:,2), 'un', false);
    actual_precomputed = CARLIN_def.seq.consites;
    verifyEqual(testCase, actual_from_bounds, expected);
    verifyEqual(testCase, actual_precomputed, expected);
end

function test_bounds_cutsites(testCase)
    check_bounds_cutsites(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_bounds_cutsites(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')));
end

function check_bounds_cutsites(testCase, CARLIN_def)    
    expected = cellfun(@(x) x(CARLIN_def.width.segment-CARLIN_def.width.cutsite+1:end), CARLIN_def.seq.segments, 'un', false);
    actual_from_bounds = arrayfun(@(s,e) CARLIN_def.seq.CARLIN(s:e), CARLIN_def.bounds.cutsites(:,1), CARLIN_def.bounds.cutsites(:,2), 'un', false);
    actual_precomputed = CARLIN_def.seq.cutsites;
    verifyEqual(testCase, actual_from_bounds, expected);
    verifyEqual(testCase, actual_precomputed, expected);
end

function test_bounds_segments(testCase)
    check_bounds_segments(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_bounds_segments(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')));
end

function check_bounds_segments(testCase, CARLIN_def)
    actual_from_bounds = arrayfun(@(s,e) CARLIN_def.seq.CARLIN(s:e), CARLIN_def.bounds.segments(:,1), CARLIN_def.bounds.segments(:,2), 'un', false);
    verifyEqual(testCase, actual_from_bounds, CARLIN_def.seq.segments);
end

function test_bounds_pams(testCase)
    check_bounds_pams(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_bounds_pams(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')));
end

function check_bounds_pams(testCase, CARLIN_def)
    actual_from_bounds = arrayfun(@(s,e) CARLIN_def.seq.CARLIN(s:e), CARLIN_def.bounds.pams(:,1), CARLIN_def.bounds.pams(:,2), 'un', false);
    verifyEqual(testCase, actual_from_bounds, CARLIN_def.seq.pam);
end

function test_bounds_postfix(testCase)
    check_bounds_postfix(testCase, CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN')));
    check_bounds_postfix(testCase, CARLIN_amplicon(parse_amplicon_file('TigreCARLIN')));
end

function check_bounds_postfix(testCase, CARLIN_def)
    verifyEqual(testCase, CARLIN_def.seq.CARLIN(CARLIN_def.bounds.postfix(1):CARLIN_def.bounds.postfix(2)), CARLIN_def.seq.postfix);
end