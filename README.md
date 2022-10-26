 TOM8 Template-O-Matic Toolbox 
 ==========================================================================

 http://irc.cchmc.org/software/pedbrain.php  
 https://neuro-jena.github.io/software.html#tom
 
 Christian Gaser christian.gaser@uni-jena.de
 Structural Brain Mapping Group, Department of Psychiatry 
 University of Jena, Germany

 Marko Wilke 
 Department of Pediatric Neurology and Developmental Medicine 
 University of Tuebingen, Germany

 Scott Holland, Mekibib Altaye  
 Imaging Research Center   
 Cincinnati Children's Hospital Medical Center, USA


Description
--------------------------------------------------------------------------

The TOM toolbox takes a radically new approach towards providing reference 
data, based on imaging data from the NIH study of normal brain development  
(http://www.bic.mni.mcgill.ca/nihpd/info/). Using the general linear model,  
we statistically isolate the influence of external variables of interest on  
brain structure, allowing us to generate high-quality matched templates for  
any given group of subjects. The toolbox offers two options:

 (1) to create pediatric templates (T1) and tissue maps (GM, WM, and CSF)  
     based on the objective 1 NIH data (n = 394), in the age range of 5-18 
     years, or
 (2) to assess a new reference population with regard to your variables of  
     interest.

Of note, this approach is generally applicable and in no way restricted to  
analyzing pediatric imaging data: for example, if you aim at investigating  
the effects of aging in elderly subjects, the toolbox will also allow you  
to create more appropriate reference (if your group is large enough to  
isolate such effects).
This toolbox is the result of a joint effort by the Department of Pediatric  
Neurology and Developmental Medicine (Marko Wilke, Tuebingen, Germany), the  
Imaging Research Center (Scott Holland and Mekibib Altaye, Cincinnati, OH,
USA), and the Structural Brain Mapping Group (Christian Gaser, Jena, Germany).
The rationale, approach and further details are available in a recent 
NeuroImage paper (Wilke et al. 2008).
We make the software available free of charge on the Imaging Research Center's
website at http://irc.cchmc.org/software/pedbrain.php.  
We will ask you to complete a registration form before downloading.

Template creation method
--------------------------------------------------------------------------

Two general approaches seem feasible to construct appropriate reference data.
First, the average age, gender, etc. is calculated based on the supplied input
information (i.e., the demographic variables of the sample under study), and a
fitting average template is created accordingly. Here, we term this the 
average approach. Alternatively, the input sample could be completely matched 
such that one reference tissue map is generated for each input subject, and 
these matched reference maps would only be averaged at the end. We term this 
the matched pairs approach.
For all template files we also save the according mat-file. Although, this is
non-standard for nifti-images, it will offer backwards compatibility to older 
SPM versions. SPM5/8 (and all other nifti-based software) will simply ignore the
mat-file.

Order of polynomial regression
--------------------------------------------------------------------------

Age can be modeled as polynomial regression with up to third order terms.
The simplest model is a linear regression (which is not recommended).
Either third or second order regression is appropriate for modeling aging
effects. The different age terms will be orthogonalized with regard to its
preceeding column.
