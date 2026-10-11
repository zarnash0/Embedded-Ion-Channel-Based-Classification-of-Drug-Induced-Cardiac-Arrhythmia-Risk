%% LINT_ALL - static check of every project .m file before a long run
files = dir(fullfile(fileparts(mfilename('fullpath')),'*.m'));
nErr = 0;
for k = 1:numel(files)
    f = files(k).name;
    if strcmp(f,'lint_all.m'), continue; end
    msgs = checkcode(f, '-struct');
    hard = msgs(arrayfun(@(m) ~isempty(strfind(lower(m.message),'parse error')) || ...
                              ~isempty(strfind(lower(m.message),'invalid')) || ...
                              ~isempty(strfind(lower(m.message),'unbalanced')) || ...
                              ~isempty(strfind(lower(m.message),'not enough')) || ...
                              ~isempty(strfind(lower(m.message),'undefined')), msgs));
    if ~isempty(hard)
        nErr = nErr + numel(hard);
        fprintf('\n### %s\n', f);
        for q = 1:numel(hard)
            fprintf('   line %4d: %s\n', hard(q).line, hard(q).message);
        end
    end
end
if nErr == 0
    fprintf('LINT CLEAN - no parse/undefined errors in any file.\n');
else
    fprintf('\n%d hard issue(s) found.\n', nErr);
end
