# Attenuated evolution of mammals through the Cenozoic

Data and novel code for Goswami et al. 2022. Attenuated evolution of mammals through the Cenozoic. Science, 378: 377-383.

## Data
- `shape.data.322.csv`: Procrustes coordinates for all specimens
- `csize.322.csv`: Centroid size for all specimens
- `full_species.data.csv`: full data on species names, clade, status, and ecological and life history traits used in analyses
- `lm_ids.csv`: list and element association of landmark and semilandmarks used in analyses
- /trees: Provides all trees used in analyses
- 3D scans available at Phenome10K.org or Morphosource, as indicated in Data S1, unless prohibited by the holding institution

## Novel code
### Trees
- Code for grafting trees is provided at: https://github.com/rnfelice/tree_graftR

### BayesTraits
- Code for generating fan plots from BayesTraits analysis is available at https://github.com/rnfelice/Dinosaur_Skulls/tree/master/Scripts). 
- `BT_results_clades`: code for generating plots comparing tip rates across clades from BT output, including generating summary and myrateslist files from BT output.

## Interactive Figures
Interactive morphospace plots corresponding to Fig 1B and Fig S2. These can be opened in a web browser. Generated with Plotly. 
