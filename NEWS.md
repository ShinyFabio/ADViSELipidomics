# ADViSELipidomics 1.3.0
* Graphical improvements
* ADViSELipidomics now can perform DA with more than two design variables uploading a txt file with the contrasts. For 1 or 2 design variables it is still possible to easily write the contrasts inside the shinyapp as before.
* Removed two package dependencies (readxl and stringi)
* Added other tooltips
* Added a docker version of ADViSELipidomics


# ADViSELipidomics 1.2.0

* Added check for the LipidSearch version.
* ADViSELipidomics now supports files from LipidSearch 5.0.
* Now it's possible to customize the height of the Lipid class proportion plot.
* Corrected a bug when the user tries to use a lipid in InternalStandardLipidIon that doesn't match the Class column (e.g. a PS(18:1D7_15:0)-H for the PG lipid class).
* Added warning messages if there are missing files in the loading process. The user will be informed which file is missing and if the target file is correctly filled (e.g. if it is used another separator instead of the semi-colon).
* Corrected a bug that forced the user to restart ADViSELipidomics if a wrong calibration file is loaded (both deuterated and nonlabeled).
* Improved the graphical interface


# ADViSELipidomics 1.1.0

* Heatmap will not cause anymore the crash of ADViSELipidomics if there are many missing values. Instead, an error pops out.
* Added more filters in the heatmap section.
* Improved visualization for many plots.


# ADViSELipidomics 1.0.0

* First release of ADViSELipidomics
