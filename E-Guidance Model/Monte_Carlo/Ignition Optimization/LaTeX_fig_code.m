% Function to generate LaTeX figure code of figures in directory
clear
close all

% inputs
output_file = 'fig_code.txt';
filetype = '.pdf';
placement = 'H';


% check if file exists. Error and exit if it does.
fid = fopen(output_file);
if fid ~= -1
    fprintf('Error: File already exists. Please specify new file\n')
    fclose(fid);
    return
else
    fclose('all')
    fid = fopen(output_file,'w');
end

% Get file names
listing = dir(strcat('*',filetype));
count = size(listing);
count = count(1);



for k = 1:count
    figname = listing(k).name(1:end-length(filetype));
fprintf(fid,'\\begin{figure}[%s]\n',placement)
fprintf(fid,'    \\centering\n')
fprintf(fid,'    \\begin{minipage}{4.5 in}\n')
fprintf(fid,'        \\includegraphics[width=\\linewidth]{Figures/%s}\n',strcat(figname,filetype))
fprintf(fid,'        \\caption{INSERT FIGURE CAPTION \\label{fig:%s} }\n',figname)
fprintf(fid,'    \\end{minipage}\n')
fprintf(fid,'\\end{figure}\n\n\n\n')

end

fclose('all')