# ADViSELipidomics 
ADViSELipidomics is a novel Shiny app for the preprocessing, analysis, and visualization of lipidomics data. It copes with the outputs from LipidSearch and LIQUID for lipid identification and quantification, and with data available from the Metabolomics Workbench. ADViSELipidomics extracts information by parsing lipid species (using LIPID MAPS classification) and, together with information available on the samples, allows performing several exploratory and statistical analyses. In the presence of internal lipid standards in the experiment, ADViSELipidomics can normalize the data matrix, providing absolute values of concentration per lipid and sample. Moreover, it allows the identification of differentially abundant lipids in simple and complex experimental designs, dealing with batch effect correction.



# 1. Install
ADViSELipidomics is a stand-alone Shiny application developed in RStudio IDE (RStudio > 1.4) and implemented using the R language (R > 4.0), available at the following GitHub page: https://github.com/ShinyFabio/ADViSELipidomics. ADViSELipidomics is multi-platform. We tested its functionalities on the main operating systems: Windows 10, Windows 11, macOS 12, Ubuntu 18, Ubuntu 20. 
The user must first install R (https://www.r-project.org) and R studio (https://www.rstudio.com), if not yet available. Then, before installing ADViSELipidomics, the user might need to perform few supplementary steps that depend on the operating systems:


ADViSELipidomics is a stand-alone application implemented using the R
language (R > 4.0) and the Shiny libraries. It can be installed as any
other R package on several operating systems (Windows, macOS and Linux).
Before installing the package you have to perform few supplementary
steps based on your operating systems:

-   **Windows** Install Rtools, a collection of tools necessary for building R packages in Windows, and it is available at the following link: <https://cran.r-project.org/bin/windows/Rtools>



-   **MacOS**  The following code should be written in the console:

``` r
brew install imagemagick@6
brew install cairo
```

-   **Ubuntu**  The following code should be written in the console:


If you are on Ubuntu run the following codes in the console:

``` r
sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
sudo apt-get install libcairo2-dev
sudo apt-get install libxt-dev
sudo apt install libmagick++-dev
sudo apt-get install libc6
sudo apt-get install libnlopt-dev
```

Then, for all the operating systems, ADViSELipidomics can be installed typing the following code in the RStudio console:

``` r
if(!require("devtools")){
  install.packages("devtools")
}
library(devtools)
install_github("ShinyFabio/ADViSELipidomics")
```

We kindly suggest updating all the R packages requested during the installation process of ADViSELipidomics Shiny application.
Finally to execute ADViSELipidomics the user can type the following code in the RStudio console:

``` r
library("ADViSELipidomics")
run_ADViSELipidomics()
```


# 2. Guide

ADViSELipidomics has a graphical user interface (GUI) implemented using the shiny and golem R packages. It has five main sections: Home, Data Import & Preprocessing, SumExp Visualization, Exploratory Analysis, and Statistical Analysis. Each section is accessible from a sidemenu on the left.


## Home section
Home section contains general information about ADViSELipidomics like the citation, the link to the github page and the link to this manual. From the "Start!" button it is possible to go to the following section where the user can upload the lipidomic data.

## Data Import & Preprocessing
This section allows users to import and process lipidomics data from various sources. 
When the user opens this section for the first time after launch, a message box appears and asks you to write your name and your company. These informations will be stored in the final output of ADViSELipidomics. By default, if you click on "Run", the following informations will be stored: User: "Name", Company: "Company". Next, the user can choose between different types of data i.e. LipidSearch, LIQUID, Excel files, Summarized Experiment, and Metabolomics Workbench. Moreover, it is also possible to select between experiments with or without internal standards (this option is available only for LipidSearch import). Based on your choice, the workflow can be different and requires different files. Here we describe all the required files for each typology.


* **LipidSearch without Internal Standard lipid: **
  * Target file (.xlsx)
  * Internal Reference file (.xlxs)
  * Output coming from LipidSearch related to your samples(.txt)
* **LipidSearch with Internal Standard lipid: **
  * Target file (.xlsx)
  * Internal Reference file (.xlxs)
  * Output coming from LipidSearch related to your samples(.txt)
  * Calibration File for deuterated (.xlsx)
  * Calibration File for nonlabeled (.xlsx)
  * Output coming from LipidSearch related to the internal standard (.txt)
* **LIQUID**
  * Target file (.xlsx)
  * Internal Reference file (.xlxs)
  * Output coming from LIQUID related to your samples(.tsv)
* **User’s Excel File**
  * Target file (.xlsx)
  * Data Matrix File (.xlsx)
* **SummarizedExperiment**
  * SummarizedExperiment object (.rds)

For **Metabolomics Workbench** you don't need to import anything, just choose the  Metabolomics Workbench ID study.


A part of the output from LipidSearch and LIQUID, ADViSELipidomics requires that the Excel files have a given structure with some mandatory columns. Here we provide a guide to the creation of these Excel files.

### Data files (LipidSearch or LIQUID)

The output of LipidSearch and LIQUID are some text files containing information on chromatographic peak area or peak intensity per lipid. If your data come from **LipidSearch** you should have a deuterated file and a non-labeled file for each sample (or replicate). The extension of these files should be .txt. If your data come from **LIQUID** you can have a positive and a negative file, with .tsv extension. In any case, put your data file in a folder and rename each file with your sample id in a proper way. 

**Example:**
your sample it's called "AF-1CM" and you have two technical replicates. Then, depending on the output software, the name of the data files should be:
* for LipidSearch: "AF-1C-M_deuterated_1.txt", "AF-1C-M_nonlabeled_1.txt" and "AF-1C-M_deuterated_2.txt" and "AF-1C-M_nonlabeled_2.txt"
* for LIQUID: "AF-1C-M_positive_1", "AF-1C-M_negative_1" and "AF-1C-M_positive_2", "AF-1C-M_negative_2"

The last two characters (e.g. "_1") refers to the technical replicate. If you don't have technical replicates just remove these two last characters (so for example in the case of LipidSearch you should have "AF-1C-M_deuterated.txt" and "AF-1C-M_nonlabeled.txt").




**NOTE**
In the choice of your sample name, it's better to avoid special characters and **DO NOT use "_"**. This character is used by ADViSELipidomics to split the file name in three parts: the sample name, the type of file (deuterated/nonlabeled or positive/negative) and the technical replicate.
Bad name: "Blood_bag_deuterated_1.txt"
Good name: "Blood-bag_deuterated_1.txt"


### Target File
The Target File is an Excel file that contains all the informations about your samples. It is the most important file since it is used for  LipidSearch import, LIQUID import, and User’s Excel File import. This file requires some mandatory columns that have to be filled with some criteria:


* **SampleID** this column contains the ID of each sample. To prevent errors, the best way to write the IDs is: "samplename_1" where in "samplename" you can write your sample name and "_1" represents the identification for the technical replicate. If your experiment doesn't have technical replicates, you can simply write "samplename". A good SampleID could be "AF-1C_1" if technical replicate are present, or "AF-1C" if not. A bad name could be "AF_1C_1" (with another underscore).


*	**File_name** this column contains the name of the data files coming from LipidSearch or LIQUID. In both cases for each sample there are two different file. In LipidSearch you have a "deuterated" and a "nonlabeled" file, while in LIQUID you have a "positive" and a "negative" file. Depending on your data type, write both file names in the corresponding cell separated by a semicolon ";".


SampleID similar to targetfile, so samplename_1 if replicates are present otherwhise jus samplename. E.g. AF-CD14-1_1
-	Norm_factor optional. If not present data won’t be normalized. Just put a number. Be careful with decimal point (if . or ,).



Here we present an example of a Target File.





ADViSELipidomics deals with the data files containing information on chromatographic peak area or peak intensity per lipid, obtained as output from LipidSearch. If your lipidomic data come from LipidSearch, you can choose to use or not internal standard
