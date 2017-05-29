https://www1.ncdc.noaa.gov/pub/data/paleo/softlib/vs-lite/readme-vslite.txt

VS-Lite Model of Tree-Ring Width 
----------------------------------------------------------------------- 
               World Data Center for Paleoclimatology, Boulder 
                                  and 
                     NOAA Paleoclimatology Program 
----------------------------------------------------------------------- 
NOTE: PLEASE CITE ORIGINAL REFERENCE WHEN USING THIS MODEL!!!!! 


NAME OF DATA SET: VS-Lite Model of Tree-Ring Width 
LAST UPDATE: 7/2013 (Option to input moisture added, rather than 
estimated from T and P)
CONTRIBUTORS: Tolwinski-Ward, S.E., M.N. Evans, M.K. Hughes, 
and K.J. Anchukaitis. 
IGBP PAGES/WDCA CONTRIBUTION SERIES NUMBER: 2010-130 

WDC PALEO CONTRIBUTION SERIES CITATION: 
Tolwinski-Ward, S.E., et al. 2010. 
VS-Lite Model of Tree-Ring Width. 
IGBP PAGES/World Data Center for Paleoclimatology 
Data Contribution Series # 2010-130. 
NOAA/NCDC Paleoclimatology Program, Boulder CO, USA. 


ORIGINAL REFERENCE: 
Tolwinski-Ward, S.E., M.N. Evans, M.K. Hughes, 
and K.J. Anchukaitis. 2011. 
An efficient forward model of the climate controls 
on interannual variation in tree-ring width. 
Climate Dynamics, 36: 11-12, 2419-2439,
doi: 10.1007/s00382-010-0945-5 

ABSTRACT: 
We present a simple, efficient, process-based forward model 
of tree-ring growth, called Vaganov-Shashkin-Lite (VS-Lite), 
that requires as inputs only latitude and monthly temperature 
and precipitation. Simulations of six bristlecone pine ring-width 
chronologies demonstrate the interpretability of model output 
as an accurate representation of the climatic controls on growth. 
Ensemble simulations by VS-Lite of two networks of North American 
ring-width chronologies correlate with observations at higher 
significance levels on average than simulations formed by 
regression of ring width on the principal components of the 
same monthly climate data. VS-Lite retains more skill outside 
of calibration intervals than does the principal components 
regression approach. It captures the dominant low- and high 
frequency spatiotemporal ring-width signals in the network 
with an inhomogeneous, multivariate relationship to climate.  
Because continuous meteorological data are most widely available 
at monthly temporal resolution, our model extends the set 
of sites at which forward-modeling studies are possible. 


ADDITIONAL REFERENCES: 
Tolwinski-Ward, S.E., K.J. Anchukaitis, and M.N. Evans. 2013.
Bayesian parameter estimation and interpretation for an intermediate
model of tree-ring width.
Climate of the Past, 9,1-13, 2013
doi: 10.5194/cp-9-1-2013

Tolwinski-Ward, S.E., M.N. Evans, M.K. Hughes, 
and K.J. Anchukaitis. 2011. 
Erratum to: An efficient forward model of the climate controls 
on interannual variation in tree-ring width. 
Climate Dynamics, 36: 11012, 2441-2445, 
doi: 10.1007/s00382-011-1062-9 

Daly, C., M. Halbleib, J.I. Smith, W.P. Gibson, M.K. Doggett, 
G.H. Taylor, J. Curtis, and P.A. Pasteris. 2008. 
Physiographically-sensitive mapping of temperature 
and precipitation across the conterminous United States. 
International Journal of Climatology, 28: 2031-2064


GEOGRAPHIC REGION: N/A 
PERIOD OF RECORD: N/A


FUNDING SOURCES: 
US National Science Foundation NSF/CMG grant 0724802, 
and NSF DMS-1204892
US National Oceanic and Atmospheric Administration grants
NOAA/CPO NA060AR4310115, NOAA NA07OAR4310060, 
and NOAA NA07OAR4310424.



DESCRIPTION: 
Vaganov-Shashkin-Lite (VS-Lite) process-based forward model 
of tree-ring growth. 

Model Revision History
 v0.1 - Original coding at monthly timestep from full daily 
        timestep model (SETW, 4/09)
 v1.0 - Changed soil moisture module to the CPC Leaky Bucket 
        model (SETW, 5/09) 
 v1.1 - No upper parametric bounds for gT, gW as in full model; 
        no density module (SETW, 9/09)
 v1.2 - Added adjustable integration window parameters (SETW, 1/10)  
 v2.0 - Minor debugging for Octave compatibility, final version 
        for publication (SETW, 10/10)
 v2.1 - Error in evapotranspiration calculation corrected (SETW, 7/11)
 v2.2 - Option for soil moisture "substepping" added (SETW, 11/11)	
 v2.3 - Add switch to allow moisture M to be given as input rather than 
        estimated from T and P; add variable input options and improved
        commenting (SETW, 7/13)

TEST RING WIDTH DATA are from the International Tree Ring Data Bank 
(http://www.ncdc.noaa.gov/paleo/treering.html)
   Site 1: 36.45N, -118.22E, 'ca530'
   Site 2: 34.17N, -117.12E, 'ca544'
   Site 3: 39.02N, -122.82E, 'ca615'
   Site 4: 39.32N, -106.08E, 'co523'

TEST CLIMATE DATA are from OSU's PRISM climate product 
(http://www.prism.oregonstate.edu/)
See Daly et al. 2008.  

If you find bugs or would like to suggest model extensions or improvements, 
please contact Suz Tolwinski-Ward, tolwinsk@ucar.edu
---------------------------------------------------------------------
VSLite package contains the following files:
1. readme-vslite.txt (this file)
2. test_vslite_v2_3.m - test code loads test data, estimates parameters, runs model, plots output.
3. VSLite_v2_3.m - main model code.
4. VSLite_v2_2_documentation.pdf - documentation for updates in version 2.2
5. vslite_testdata.mat - file of test data 
6. VSLite_correction_7_2011.pdf - documentation for correction to model code
7. possible-erratum-email.txt - regarding typographical error (fixed v2.1) in code
8. estimate_vslite_params_v2_3.m - parameter estimation code for VSLite v2.3



To test VS-lite code, type

test_vslite_v2_3

at the MATLAB command prompt. Script is set up to load test data from vslite_testdata.mat, call parameter estimation code to estimate parameters (with default options), and to call the main code of VS-Lite to run test simulations (with all default options) and finally display output. 

Typing 

help VSLite_v2_3

and

help estimate_vslite_params_v2_3

at the MATLAB command prompt will give more general instructions and list all user-options for running the VS-Lite model and the parameter estimation code.

SETW 7/9/2013
