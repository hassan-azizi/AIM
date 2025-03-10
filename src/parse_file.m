function file_data = parse_file(filename)
    % Function to parse *.csv or *.dat or *.txt files containing isotherm
    % data
    
    % Initialize the struct object
    file_data = struct();

    file_data.meta_data = dictionary;
    file_data.data = [];
    
    % Matching Pattern
    match_pat = regexpPattern('(\-\d*|\d*)\.\d*|(\-\d*|\d*)');
    
    % Opening the file
    fileID = fopen(filename, "r");
    
    %% Looping every line of file
    while ~feof(fileID)        
        current_line = fgetl(fileID);

        % Trimming whitespaces
        current_line = strip(current_line);
    
        % Ignoring empty line
        if isempty(current_line)
           continue
        end
        %% Handling meta data
        if startsWith(current_line, '#')
            tokens = split(string(current_line));
            key = tokens(1);
            value = join(tokens(2:end), ' ');
            file_data.meta_data(key) = value;
            continue
        end

        %% The loop will terminate when the first line with numerical data is find
        if startsWith(current_line, match_pat)
            break
        end
    end

    % Closing the file
    fclose(fileID);
    
    isotherm_data = readmatrix(filename);
    
    file_data.data = isotherm_data;
end