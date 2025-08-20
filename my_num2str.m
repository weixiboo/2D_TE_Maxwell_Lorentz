function latexStr = my_num2str(number)
    if number >= 1e4 || number < 1e-2
        scientificStr = num2str(number, '%.3e'); % Example: '1.234568e+05'
        latexStr = strrep(scientificStr, 'e+0', ' \times 10^{');
        latexStr = strrep(latexStr, 'e-0', ' \times 10^{-');
        
        latexStr = strrep(latexStr, 'e+', ' \times 10^{');
        latexStr = strrep(latexStr, 'e-', ' \times 10^{-');
    
        latexStr = strcat(latexStr,'}$');
        latexStr = strcat('$',latexStr);
    else
        number = round(number,3);
        latexStr = num2str(number);
    end

     if isnan(number)
        latexStr = '-';
    end

end

