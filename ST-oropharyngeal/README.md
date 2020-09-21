# Oropharyngeal 

## Data description
**Data acquisition**
All patients provided informed consent for the collection of human specimens and data.
This was approved by the St Vincent’s Hospital Research Office (2019/PID04335) in accordance with the National Health and Medical Research Council’s National Statement of Ethical Conduct in Human Research.  Patients undergoing surgical resection for a locally advanced oropharyngeal cancer were recruited to the study.  After surgical removal, the anatomical pathologist dissected a sample of both the primary and nodal metastasis.
Samples were tumour banked in accordance with our ethically approved protocol.

**Sample storage**
Within 30 minutes of collection, tumour samples were tumour banked.  Samples were cut into 1mm x 1mm chunks with a scalpel blade.
For Visium, a tissue chunk was snap frozen in OCT. After freezing, samples were moved to liquid nitrogen for long term storage.

**Visium Spatial Gene Expression**
Frozen tissue samples were processed using the Visium Spatial Gene Expression slide and reagent kit (10X Genomics, US) following the manufacturer’s instruction. 
Briefly, 10 μm sections were placed into the capture areas of the Visium slide.
Tissue morphology was assessed with H&E staining and imaging using a Leica DM6000 microscope equipped with a 20x lens (Leica, DE).
The imaged sections were then permeabilized for 12 minutes using the supplied reagents.
The permeabilization condition was previously optimised using the Visium Spatial Tissue Optimisation slide and reagent kit (10X Genomics, US).
After permeabilization, cDNA libraries were prepared, checked for quality and sequenced on a NovaSeq 6000 platform (Illumina, US).
Around 300 million pair-ended reads were obtained for each tissue section. Read 1, i7 index and Read 2 were sequenced with 28, 8 and 98 cycles respectively.

## Data availability
Pending to get approval to post to GEO.

## Code
Scripts 1-_australia_oroph_processing.Rmd, 2-australia_oroph_biological.Rmd, and 3-australia_oroph_deconv.Rmd are in charge of preprocessing and mapping 
the TICA immune cell states to the tissue slices. 4-australia_oro_srtatification.Rmd  and 5-australia_oro_plots.Rmd, in turn, are in charge of making the plots for Figure 5 and Supplementary Figures 8-12.
All code was run with R 3.6.3.
