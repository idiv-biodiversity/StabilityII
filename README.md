---
output: 
  html_document: 
    theme: flatly
---
## Code to support analyses in:    
Authors: Craven, D., N. Eisenhauer, W.D. Pearse, Y. Hautier, F.Isbell,C. Roscher, M. Bahn, C.
Beierkuhnlein, G. Bönisch, N. Buchmann, C. Byun, J.A. Catford, B.E.L. Cerabolini, J.H.C.
Cornelissen, J.M. Craine, E. De Luca, A. Ebeling, J.N. Griffin, A. Hector, J. Hines, A.
Jentsch, J.Kattge, J. Kreyling, V. Lanta, N. Lemoine, S.T. Meyer, V. Minden, V. Onipchenko,
H.W. Polley, P.B. Reich, J. van Ruijven, B. Schamp, M.D. Smith, N.A. Soudzilovskaia, D.
Tilman, A. Weigelt, B. Wilsey, and P. Manning. 2018. **Multiple facets of biodiversity
drive the diversity-stability relationship**. _Nature Ecology & Evolution_ In revision.

**Data is currently available upon request. Please send an e-mail to:** dylan.craven@aya.yale.edu


## Table of Contents  

1. **LME_BiodivFacets_Stability.R** : linear mixed-effects models testing variation in ecosystem stability for each facet of biodiversity, i.e. species richness,species asynchrony, phylogenetic diversity, fast-slow functional diversity, and community-weighted means of fast-slow traits.

2. **Pairwise_Correlation.R**: Pairwise correlations among facets of biodiversity using multi-level meta-analytical regression. 

3. **TS_basicSEM.R**: Basic piecewise SEM for each combination of phylogenetic and functional diversity metric using all data (4 models).

4. **TS_basicSEM_longterm.R**: Basic piecewise SEM for _only_ long-term studies (> 4 years) using each combination of phylogenetic and functional diversity metrics (4 models).  

5. **TS_basicSEM_indtrts.R**: Basic piecewise SEM for each individual functional trait, i.e. SLA, LDMC, leaf N, leaf P, using each combination of phylogenetic and functional diversity metrics (16 models: 4 combinations of functional and phylogenetic diversity metrix x 4 individual functional traits). 

6. **TS_extendedSEM.R**: Extended piecewise SEM, i.e. including mean and SD of aboveground biomass production, for each combination of phylogenetic and functional diversity metrics (4 models).  

7. **Calculate_DetrendedStability.R**:  Calculate detrended ecosystem stability using linear regression for each plot.  

8. **LME_BiodivFacets_DetrendedStability.R** : linear mixed-effects models testing variation in **detrended ecosystem stability** for each facet of biodiversity, i.e. species richness,species asynchrony, phylogenetic diversity, fast-slow functional diversity, and community-weighted means of fast-slow traits. 

9. **TS_basicSEM_DetrendedStability.R**: Basic piecewise SEM for each combination of phylogenetic and functional diversity metric using all data (4 models) and **detrended ecosystem stability**.

10. **TS_extendedSEM_DetrendedStability.R**: Extended piecewise SEM, i.e. including mean and SD of aboveground biomass production, for each combination of phylogenetic and functional diversity metrics (4 models) and **detrended ecosystem stability**.  