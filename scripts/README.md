## Documentation Scripts
Author: Kelsey Keith ([@kelseykeith]())

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Purpose](#purpose)
- [Usage](#usage)
- [Results](#results)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

### Purpose

In the Pediatric Molecular Open Targets high-level documentation, the main README of this repository, we have a table with the number of samples per cohort for the Datasets section so that people are aware of the general number of samples they can access. The script `biospecimen_count.Rmd` calculates and outputs this table and well as a table giving the sample counts per disease and dataset for reference.

Data can be accessed through the `OpenPedCan-analysis` repository <https://github.com/PediatricOpenTargets/OpenPedCan-analysis>, although for privacy reasons not all data, especially raw data is accessible publically.

### Using the `OpenPedCan-analysis` Submodule

In order to get the data in the `OpenPedCan-analysis` submodule you have to first initialize the submodule, then use the `data-download.sh` script within the submodule to fetch the data using the following commands:

```
# initialize the submodule
# Note: you should see a message that the submodule has been initialized followed by the standard git progress messages
git submodule update --init --progress

# change directory into the submodule and download the data
cd OpenPedCan-analysis
bash download-data.sh
``` 

### Usage

You can either open `biospecimen_count.Rmd` in Rstudio and knit it or execute it on the command line with the following command:

```bash
Rscript -e "rmarkdown::render('biospecimen_count.Rmd', clean = TRUE)"
```

### Results

The script outputs three files:

- `./biospecimen_dataset_counts.csv`: a plain text version of the table pictured in the main repository README giving a count of the DNA and RNA biospecimens in each Dataset
- `../figures/biospecimen_dataset_counts.png`: a nicely formated image of the table that is inserted in the main README giving a count of the DNA and RNA biospecimens in each Dataset
`./biospecimen_dataset_disease_counts.csv`: a plain text table giving a count of the DNA and RNA biospecimens in each Disease in each Dataset for reference

<br><br>