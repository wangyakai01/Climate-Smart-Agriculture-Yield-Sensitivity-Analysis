
function coefTableWithStars=extract_fixed_effect_Res(T)
numModels = height(T);
allTerms = {};

% Step 1: collect all unique fixed-effect term names
for i = 1:numModels
    mdl = T.lme{i};
    coefs = mdl.Coefficients;
    allTerms = union(allTerms, coefs.Name);  % gather all fixed effect terms
end

numTerms = length(allTerms);
coefStrings = strings(numModels, numTerms);  % for final display

% Step 2: fill in estimates + significance markers
for i = 1:numModels
    mdl = T.lme{i};
    coefs = mdl.Coefficients;
    
    for j = 1:height(coefs)
        term = coefs.Name{j};
        est = coefs.Estimate(j);
        pval = coefs.pValue(j);
        
        % Significance symbols
        if pval <= 0.01
            mark = '**';
        elseif pval <= 0.05
            mark = '*';
        else
            mark = '';
        end
        
        idx = find(strcmp(allTerms, term));
        coefStrings(i, idx) = sprintf('%.3f%s', est, mark);
    end
end

% Step 3: create table
coefTableWithStars = array2table(coefStrings, 'VariableNames', matlab.lang.makeValidName(allTerms));
coefTableWithStars.Crop = T.crop;
coefTableWithStars = movevars(coefTableWithStars, 'Crop', 'Before', 1);

% Step 4: show the table
disp(coefTableWithStars);
end