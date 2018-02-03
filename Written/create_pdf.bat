cd C:\Users\stroll\Google Drive\Classes\SDSU\Thesis\Written
latex thesis
bibtex thesis
latex thesis
latex thesis
dvips -o thesis.ps thesis
epstopdf --outfile thesis.pdf  thesis.ps