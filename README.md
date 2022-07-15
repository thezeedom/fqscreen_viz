 R Shiny app for visualizing FastQ Screen data.

- app.R
- fqscreen_function_MASTER.R - Functions utilized in app.

## Input Data  
This app requires the text output from <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/" target="_blank">FastQ Screen</a>.

The app parses filenames to determine the unique samples and the type of read (either single-end or paired-end).

For paired-end, it is required for paired sample file names to include reference to R1 and R2, respectively.

## Output

Files are melted to data.frames containing percent mapping data and number mapping data. Percent mapping data is graphed per sample.
