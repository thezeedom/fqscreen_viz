 R Shiny app for visualizing FastQ Screen data.

- app.R - Main application. Recommeded to run in browser.
- fqscreen_function_MASTER.R - Functions utilized in app.

## Input Data  
This app requires the text output from <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/" target="_blank">FastQ Screen</a>.

To include multiple files, please run the app in browser. 

The app parses filenames to determine the unique samples and the type of read (either single-end or paired-end). To this end it is assumed that filenames include a reference to either R1, R2 or neither. 

*Example filenames:*  
1_S9_R1_screen.txt  
129_S19_screen.txt  
WTIFNB_1_S7_R1_screen.txt  
Palb-Tram-TTR-20R_S49_screen.txt  
Ulix-23L_S45_screen.txt

## Output

Files are melted to data.frames containing percent mapping data and number mapping data. Percent mapping data is used for plotting.