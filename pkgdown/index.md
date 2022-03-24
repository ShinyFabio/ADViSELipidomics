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



# 2. Input Data

ADViSELipidomics allows the user to import files concerning different types of data:
* LipidSearch or LIQUID. ADViSELipidomics deals with the data files containing information on chromatographic peak area or peak intensity per lipid, obtained as output from external software for identifying and quantifying lipids (i.e., ADViSELipidomics currently supports the output formats from LipidSearch or LIQUID). Moreover, it requires the Target File with details on samples (such as treatments or biological replicates), the Internal Reference File with bounds for the filtering step in the following modules. ADViSELipidomics shows a quality plot based on the sum of chromatographic peak area per sample (or replicate). In the case of LipidSearch output associated with internal lipid standards, ADViSELipidomics also requires all the Calibration Files for the construction of the calibration curves.
* Metabolomics Workbench. ADViSELipidomics can download in real-time suitable selected lipidomic experiments from the online repository;
* Excel. The user can upload two Excel files, the data matrix and the Target File;
* SummarizedExperiment. The user can upload a SummarizedExperiment R object (SE), with several types of information (data matrix, information on lipids, information on samples, metadata if available).

Hence, as can be seen, ADViSELipidomics requires different files that may change beetween the different types of data. To sum up, here is a list with all the required files for each data type:


* **LipidSearch without Internal Standard lipid:** 
  * Target file (.xlsx)
  * Internal Reference file (.xlxs)
  * Data files coming from LipidSearch related to your samples (.txt)
* **LipidSearch with Internal Standard lipid:** 
  * Target file (.xlsx)
  * Internal Reference file (.xlxs)
  * Data files coming from LipidSearch related to your samples (.txt)
  * Calibration File for deuterated (.xlsx)
  * Calibration File for nonlabeled (.xlsx)
  * Concentration files coming from LipidSearch related to the internal standard (.txt)
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

Before run ADViSELipidomics make sure that you have all the required files and that they are draw up properly. Apart of the output files from LipidSearch and LIQUID, ADViSELipidomics requires that the Excel files have a given structure with some mandatory columns. Here we provide a guide to the creation of these Excel files. 



## 2.1 Data files (LipidSearch or LIQUID)

The output of LipidSearch and LIQUID are some text files containing information on chromatographic peak area or peak intensity per lipid. If your data come from **LipidSearch** you should have a deuterated file and a non-labeled file for each sample (or replicate). The extension of these files should be .txt. If your data come from **LIQUID** you can have a positive and a negative file, with .tsv extension. In any case, put your data file in a folder and rename each file with your sample id in a proper way. 

**Example:** 
your sample it's called "AF-1CM" and you have two technical replicates. Then, depending on the output software, the name of the data files should be:
* for LipidSearch: "AF-1C-M_deuterated_1.txt", "AF-1C-M_nonlabeled_1.txt" and "AF-1C-M_deuterated_2.txt" and "AF-1C-M_nonlabeled_2.txt"
* for LIQUID: "AF-1C-M_positive_1", "AF-1C-M_negative_1" and "AF-1C-M_positive_2", "AF-1C-M_negative_2"

The last two characters (e.g. "_1") refers to the technical replicate. If you don't have technical replicates just remove these two last characters (so for example in the case of LipidSearch you should have "AF-1C-M_deuterated.txt" and "AF-1C-M_nonlabeled.txt").



**NOTE**
In the choice of your sample name, it's better to avoid special characters and **DO NOT use underscores (_)**. This character is used by ADViSELipidomics to split the file name in three parts: the sample name, the type of file (deuterated/nonlabeled or positive/negative) and the technical replicate as shown in the following picture:

![datafile_names](https://user-images.githubusercontent.com/78078351/159508536-ee7bc3c2-59c2-4c11-bb32-782c7eb76696.png)

Bad name: "Blood_bag_deuterated_1.txt"
Good name: "Blood-bag_deuterated_1.txt"


## 2.2 Target File (LipidSearch, LIQUID, and User’s Excel File)
The Target File is an Excel file that contains all the informations about your samples. It is the most important file since it is used for  LipidSearch import, LIQUID import, and User’s Excel File import. This file requires some mandatory columns that have to be filled with some criteria:


* **SampleID (LipidSearch, LIQUID and User’s Excel)** this column contains the ID of each sample. To prevent errors, the best way to write the IDs is: "samplename_1" where in "samplename" you can write your sample name and "_1" represents the identification for the technical replicate. If your experiment doesn't have technical replicates, you can simply write "samplename". A good SampleID could be "AF-1C_1" if technical replicate are present, or "AF-1C" if not. A bad name could be "AF_1C_1" (with another underscore).
*	**File_name (LipidSearch, LIQUID)** this column contains the name of the data files coming from LipidSearch or LIQUID. In both cases for each sample there are two different file. In LipidSearch you have a "deuterated" and a "nonlabeled" file, while in LIQUID you have a "positive" and a "negative" file. Depending on your data type, write both file names in the corresponding cell separated by a **semicolon ";"**.
* **Norm_factor (LipidSearch, LIQUID - optional)** If you need to normalize your data by a normalization factor, you can add this column and write a number (Be careful with decimal point) for each sample. If it's not present, data won’t be normalized.

In the picture below there is a Target File example, enlightned in yellow the mandatory columns, and in green the optional column. You can fill the Target File with any other informative column, just try to avoid special characters like \^$.?*|+()[]{} and whitespace.You can use - or _ instead of whitespace.

![Screenshot (197)](https://user-images.githubusercontent.com/78078351/159510916-f0fac7fa-fa98-4bdb-a88a-505c0f6a6098.png)

The example target file can be downloaded from here:

[Targetfile_Lipidomics.xlsx](https://github.com/ShinyFabio/ADViSELipidomics/files/8325363/Targetfile_Lipidomics_new.xlsx)


##  2.3 Internal Reference File (LipidSearch, LIQUID)
In LipidSearch and LIQUID option, ADViSELipidomics requires also another Excel file here called Internal Reference File which contains the list of the Internal Standard lipids defined per class and adduct, upper/lower bounds for the number of carbon atoms, upper/lower bounds for the number of double bonds, nominal standard concentration, and upper/lower bounds for the concentration linearity in the calibration curves. This file has many mandatory columns that depends both on the external software (LipidSearch, LIQUID) and the presence of internal standard (only for LipidSearch).

* **LipidSearch** 
  - **Class** lipid class of interest according to the nomenclature of LipidSearch (e.g. *DG* )
  - **Ion** the ion of interested written according to the nomenclature of LipidSearch(e.g. *M-H* )
  - **MinRt** minimum retention time of the class (a number)
  - **MaxRt** maximum retention time of the class (a number)
  - **InternalStandardLipidIon** lipid of .........????? (e.g. *Cer(d18:1_17:0)-H* )
  - **MinLinearity** minimum value for the range of linearity in the calibration curves. (a number) ONLY IF YOU USE INTERNAL STANDARD
  - **MaxLinearity** maximum value for the range of linearity in the calibration curves. (a number) ONLY IF YOU USE INTERNAL STANDARD

The picture below shows an Internal Reference File example in the case of LipidSearch and the presence of Internal Standards. In yellow the mandatory columns, and in green the columns needed only in the presence of Internal Standards.

![Screenshot (199)](https://user-images.githubusercontent.com/78078351/159517154-aca32f9c-1c0f-4299-9937-b49fa1b3b028.png)

The Internal Reference File example for the LipidSearch with Internal Standards option can be downloaded from here:

[Internal_Reference_file_LipidSearch_withIS.xlsx](https://github.com/ShinyFabio/ADViSELipidomics/files/8325722/Internal_Reference_file_LipidSearch_withIS.xlsx)

* **LIQUID**
  - **Class** lipid class of interest according to the nomenclature of LIQUID (e.g. *DG* )
  - **Adduct** the ion of interested written according to the nomenclature of LIQUID (e.g. *[M-H]+* )
  - **MinRt** minimum retention time of the class (a number)
  - **MaxRt** maximum retention time of the class (a number)

The picture below shows an Internal Reference File example in the case of LIQUID. In yellow the mandatory columns.

![Screenshot (201)](https://user-images.githubusercontent.com/78078351/159518027-cb2a3eea-923f-49b7-8633-0502d0935032.png)

The Internal Reference File example for the LIQUID option can be downloaded from here:

[Targetfile_Lipidomics_LIQUID.xlsx](https://github.com/ShinyFabio/ADViSELipidomics/files/8325726/Targetfile_Lipidomics_LIQUID.xlsx)


## 2.4 Calibration Files (LipidSearch with Internal Standards)
In the case of LipidSearch, if you have Internal Standard, you can choose to use or not them. In this case, you need to upload also some Calibration Files which are two Excel Files and the data files coming from LipidSearch (here called concentration files). The concentration files are the same .txt files described in chapter 2.1. Please, refer to that chapter if you need more informations about how to rename the files. Be sure that all the concentration files are inside a folder and they aren't mixed with the data files of the chapter 2.1. 
Next, ADViSELipidomics, requires two Calibration Excel files, one for the Nonlabeled and the other for the Deuterated. They share the same structure:

* **Concentration (ng/mL)** the concentration of the standard
* **Class** the lipid classes of interest separated by a comma **,** (e.g. *PG,PS,PI,PE,SM,PC,TG,DG* )
* **Name** the name of the data files coming from LipidSearch. They have to match perfectly with the file names. If you have technical replicates, separate them by a comma **,** (for example in the deuterated: *ISMix_5ugmL_deuterated_1,ISMix_5ugmL_deuterated_2,ISMix_5ugmL_deuterated_3*)

The picture below shows an example of a Calibration Excel file for the deuterated. 


![Screenshot (203)](https://user-images.githubusercontent.com/78078351/159525426-c776cb28-fa19-420b-9569-2143775edca5.png)

A toy example for the Calibration Deuterated and Calibration Nonlabeled Excel files can be downloaded here:

[Calibration_Deuterated.xlsx](https://github.com/ShinyFabio/ADViSELipidomics/files/8326063/Calibration_Deuterated.xlsx)

[Calibration_NonLabeled.xlsx](https://github.com/ShinyFabio/ADViSELipidomics/files/8326075/Calibration_NonLabeled.xlsx)


## 2.5 User’s Excel File

If you already have a  matrix file containing the abudance for each lipid, you need just two Excel files: the Target File and the Data Matrix File. Here the Target File has only one mandatory column, the SampleID. The Data Matrix (.xlsx file) must have the list of the lipids in the first column, that must be called *"Lipids"*, and then the samples (or replicates) in the following columns, with the column names that is the same of the *SampleID* of the Target File. It's not necessary that the matrix is full (i.e. without missing values) since after uploaded, it's possible to filter and impute NAs. The picture below shows an example of the Data Matrix.

![Screenshot (204)](https://user-images.githubusercontent.com/78078351/159945580-ea7466b3-fb11-4e19-87a3-e3eb6a64eeb3.png)

**NOTE:** the column names in the data matrix must follow the same rules of the SampleID for every Target File. Check chapter 2.2.

## 2.6 SummarizedExperiment

ADViSELipidomics allows the user to load a SummarizedExperiment (SE) object, saved as a .rds file, already prepared or previously downloaded after running ADViSELipidomics. Since the required SE object has a complex structure, we do not recommend the user to upload a SE object that wasn't download from ADViSELipidomics. The idea behind this option was that the user can save the SE object after the preprocessing steps and performs the exploratory and statistical analysis in another moment.


## 2.7 Metabolomics Workbench

In the case of Metabolomics Workbench, you don't need to import anything, because ADViSELipidomics downloads a selected Metabolomics Workbench experiment and converts it into an SE object.

# 3. Guide

ADViSELipidomics has a graphical user interface (GUI) implemented using the shiny and golem R packages. It has five main sections: Home, Data Import & Preprocessing, SumExp Visualization, Exploratory Analysis, and Statistical Analysis. Each section is accessible from a sidemenu on the left.


## 3.1 Home section
Home section contains general information about ADViSELipidomics like the citation, the link to the github page and the link to this manual. From the "Start!" button it is possible to go to the following section where the user can upload the lipidomic data.


## 3.2 Data Import & Preprocessing
This section allows to import and process lipidomics data from various sources. 
When you open this section for the first time after launch, a message box appears and asks you to write your name and your company. These informations will be stored in the final output of ADViSELipidomics. By default, if you click on "Run", the User will be *"Name"* and the Company will be *"Company"*. 
The picture below shows the Data Import & Preprocessing section (with the different parts enlightened with red rectangles). 


![workflow_withbox](https://user-images.githubusercontent.com/78078351/159950379-e3885176-0ccc-4bf1-b1eb-9f4378d468e2.png)

**Rectangle A** allows the user to choose between LipidSearch, LIQUID, Excel files, Summarized Experiment, and Metabolomics Workbench. Moreover, it is also possible to select between experiments with or without internal standards (this option is available only for LipidSearch import). **Rectangle B** shows three different steps for Importing & Filtering: importing data, storing and reading data, filtering data. **Rectangle C** shows five additional steps for Calibration: importing calibration files, storing calibration files, selection of the folder for the results, selection of calibration options, application of recovery. Finally, **rectangle D** shows two different steps for Filtering and Missing Data imputation and creating the SummarizedExperiment object. Note that the layout of the Data Import & Preprocessing section and the required files depends on the type of input data format that the user chooses. Go to chapter 2 if you need help to gather all the required files.

Since the option LipidSearch output with Internal Standard (IS) has the largest number of required file and of steps, here we provide a complete guide for this case. Anyway, this guide applies also to LipidSearch without IS and to LIQUID: in these cases the only difference is that there isn't the CALIBRATION module (chapter 3.2.2).

### 3.2.1 LipidSearch (IS) EXAMPLE - IMPORTING & FILTERING module
The first module is the IMPORTING & FILTERING module where the user can upload the Target File, the Internal Reference File and the Data file that come from LipidSearch.

![import filtering_module](https://user-images.githubusercontent.com/78078351/159956088-cb4bff48-7b78-4005-9dba-016ecca70e53.png)

* **Step 1.** The first files that you have to import are the Target File and the Internal Reference File (Rectangle A, steps 1 and 2). Next to each of them there is a button (yellow squared rectangle) that allows you to edit the Excel files. You can select only the columns that you need, filter the rows by one or more conditions and download the edited data. If you need to edit them, to apply the editing you have to enable the button next to the download button and click on the "Done" button (right top corner). Anyway an help button guide you through the editing options.
* **Step 2.** Here you choose the folder containing the data files coming from LipidSearch (only the data files related to the samples and NOT to the IS). After selected the folder, click on the "Read Data" button and ADViSELipidomics will start to read all the data files. A progress bar shows the percentage of completion. When the reading process is completed, you can perform a quality check on the area of each sample by clicking to the "Quality check" button.
* **Step 3.** Finally, here you can filter non-informative lipids based on retention time in the range, number of carbon atoms in the range, even number of carbon atoms, number of double bonds in the range, duplicated lipids. From the two sliders you can choose the range for the carbon number and the double bound number, while the other filters come from the Internal Reference File. If there are duplicated lipids (same m/z values for lipid peaks), ADViSELipidomics takes only the lipids with the maximum peak area. The "Filter Data" button starts this process. At the end, you can check the filtered data for each sample.

### 3.2.2  LipidSearch (IS) EXAMPLE - CALIBRATION module
If the previous module is successfully completed, the CALIBRATION module appears next to it. The Calibration module creates the calibration curves and the calibration matrix. It uses the Internal Lipid Standards reported in the Internal Reference file, and the correspondence between the Concentration Files and the lipid classes declared in the Calibration File. This module extracts the relationships between peaks area and concentration values for each internal lipid standard, constructing the calibration curves with a linear model and plotting them. The linear regression model can be classical or robust, with zero or non-zero intercept. Finally, the calibration matrix resumes all the points from the calibration curves. After the calibration process, ADViSELipidomics stores slope and intercept values for the recovery module.
As already stated, this module appears only if you are using LipidSearch output with Internal Standard (and you clicked on "Yes" in the radiobutton that asks *"you Do you have internal standard?"*). In this module, you need two Calibration Files (.xlsx, see chapter 2.4) and the Concentration files coming from LipidSearch related to the internal standard (.txt, see chapter 2.1). 

![CALIBRATION module](https://user-images.githubusercontent.com/78078351/159966937-26f75ea8-5ac0-405c-a92b-019f868d06c5.png)

* **Step. 1** Here you can upload the Calibration Files (.xlsx) for both Deuterated and Nonlabeled. This step is very similar to the Step 1. of the IMPORTING & FILTERING module. You can find further informations about these files in the chapter 2.4.
* **Step. 2** Select the folder containing the Concentration files coming from LipidSearch and then click on "Read the concentration files". Also this step is very similar to the step 2 of the previous module.
* **Step. 3** Here you can select the folder where saving the output from LipidSearch. ADViSELipidomics creates for you the folder structure.
* **Step. 4** In this step you can choose some calibration options and visualize the calibration plot for each standard.
* **Step. 5.** Finally, you can apply the recovery percentage on the concentration values for each lipid, considering the Internal Lipid Standards as lipid class reference. This normalization provides absolute concentration values for the lipids and the resulting concentration matrix can be seen by clicking on the "Check concentration matrix" button. Here it's possibile also to visualize the missing values(if applicable). Moreover, from the "Download LOL" button you can download a table containing all the lipids filtered because outside of the linearity range. 


### 3.2.3 LipidSearch (IS) EXAMPLE - MISSING DATA & SUMMARIZED EXPERIMENT
This is the last module of the preprocessing menu where you can filter and impute missing values (NAs), build the SummarizedExperiment object and download it. 

![MISSING DATA & SUMMARIZED EXPERIMENT](https://user-images.githubusercontent.com/78078351/159970295-85bcb1fe-9fb1-4061-aa98-905c78933330.png)

* **Step 1.** In the first step ADViSELipidomics computes the percentage of NAs, for each lipid (matrix rows) and each sample (matrix columns). Second, it allows retaining only lipids and/or samples with a percentage of missingness below thresholds chosen using the sliders. For example if you set *Max missing data percentage allowable on lipids* to 0.3 *Max missing data percentage allowable on samples* to 0.6 that means that only lipids (rows) with less than 30% of NAs and samples (columns) with less than 60% of NAs are stored. After that, by clicking on the "Check filtered NAs" button, ADViSELipidomics provides the missing data distributions and the data dimension before and after filtering NAs. 
* **Step 2.** Next you can impute the remaining NAs with different imputation methods, three Not Model Based (mean, median and knn) and one Model Based (irmi).
* **Step 3.** In the final step ADViSELipidomics build the SummarizedExperiment object and download it. By clicking on the "See the results" you will be redirected to the next menu, SumExp Visualization.

## 3.3 SumExp Visualization
Once you ended successfully the Preprocessing module, the first thing that you can do is check the just created SummarizedExperiment (SE) object. It can be done in SumExp Visualization menu. The complex structure of the SE object can be investigated by a red gear icon where you can choose what part of the SE object should be shown and to summarise the data (if you have technical replicates).

![SumExp Visualization](https://user-images.githubusercontent.com/78078351/159973518-6b7fc106-020a-4185-8595-ff611671868d.png)

In the picture above, it is shown the rowData part of the SE object containing the annotation on lipids. Each lipid in the "Lipids" contains an hyperlinks to the SwissLipids online repository to provide structural, biological, and analytic details.

## 3.4 Exploratory Analysis
The Exploratory Analysis menu includes three sub-menus: Plots, Clustering, and Dimensionality Reduction. 

### 3.4.1 Plots
This sub-menu allows the user to create different types of plots to show trend and behavior of data, exploring them from lipid and/or sample points of view. It has 4 panels: Lipids, Scatterplot, Heatmap and Quality plots.

* **Lipid plots.** It is possible to 1) represent the lipid class distribution (counts of lipids per class) with a pie chart, boxplot, and spider plot; 2) visualize the percentage proportion of lipid class for each sample using a barplot, 3) compare the lipid species abundance for each condition; 4) inspect the abundance of a lipid, selected by the user, in relationship with a feature from the target file (e.g., treatment) using boxplots;
* **Scatterplots.** It is possible to visualize the relationship between lipid abundance in two samples;
* **Heatmap.** It provides a highly customizable heatmap to show possible clusters among lipids or samples. The user can select many parameters: a) row annotation with the feature from the target file, b) column annotation with the information from lipids parsing, c) dendrograms for lipids and/or samples, d) distance function (Euclidean, maximum, Canberra), e) clustering method (complete, average, median, Ward); f) number of clusters for lipids and/or samples. The user can select an area in the overall heatmap and have a detailed zoom of the area itself, with associated information;
* **Quality plots.** It provides different typologies of plots (barplot, boxplot, density plot) to show the total amount of abundance (logarithmic scale) per sample, considering as reference a feature from the target file, to show possible unexpected behavior among samples or among replicates for the same sample.

The picture below shows an example of the Lipid class proportion plot.
![Lipid class proportion](https://user-images.githubusercontent.com/78078351/159976950-0e3220d2-dc21-40a2-a2a6-4256b30a5d14.png)



### 3.4.2 Clustering
The Clustering sub-menu allows the user to cluster the data by lipids or samples. The user can choose the number of clusters and the clustering method among the following algorithms: hierarchical clustering (using single, complete, Ward as linkage function) or partitioning clustering (k-means, PAM, Clara). If you choose a partitioning clustering, ADViSELipidomics performs first a PCA. Additional plots, such as the silhouette plot, can suggest the number of clusters to use.

![Clustering](https://user-images.githubusercontent.com/78078351/159980474-1cc6dbe4-e745-44d5-87ce-72f636454d2a.png)


### 3.4.3 Dimensionality Reduction
The Dimensionality Reduction sub-menu allows the user to choose between unsupervised (PCA) and supervised approaches (PLS-DA, sPLS-DA) to represent the data in a two or three-dimensional space. It contains three panels PCA, PLS-DA, and sPLS-DA.

* **PCA.** ADViSELipidomics computes the Principal Component Analysis (PCA), showing the results with different plots: a) 2D plot, b) biplot, c) scree plot, d) loadings plot, e) 3D plot. The user can highlight the features from the target file with different colors and select the number of components to use for the loading plots.

![PCA](https://user-images.githubusercontent.com/78078351/159980009-d4aa02bc-87a0-4bb0-a5af-8cc8f4df0171.png)

* **PLS-DA.** ADViSELipidomics computes the Partial Least Square - Discriminant Analysis,showing the results with a 2D plot and a Correlation Circle plot. The user can select the group variable and the number of components for the computation. The 2D plot can be customized from the red gear icon. Furthermore, ADViSELipidomics can perform a Cross-Validation in order to identify the best number of components. It may take a while.
* **sPLS-DA.** ADViSELipidomics can compute also the sparse version of the PLS-DA. The panel is very similar to the PLS-DA panel but since it's a sparse version, it is possible to choose the number of variables to select on each component (called "KeepX"). Also here ADViSELipidomics can perform a Cross-Validation that helps the user to choose the best number of components and the best "KeepX". 


## 3.5 Statistical Analysis

