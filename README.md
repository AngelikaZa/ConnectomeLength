# ConnectomeLength

Code accompanying:
    Zarkali et al. Dementia risk in Parkinson’s disease is associated with interhemispheric connectivity loss and determined by regional gene expression. NeuroImage Clinical (in       press)
    
<b>Usage</b>:

1) Data folder: 
Contains demographics, connectome level data for participants and module allocations derived from community Louvain algorithm (gamma=1).
- Demographics: Participant.csv
- Normalised connectome (glasser atlas): Data/Connectomes/<Participant_ID>/connectome_norm.csv 
- Streamline length-transformed connectome (glasser atlas): Data/Connectomes/<Participant_ID>/new_connectome.csv
- Module allocations: These are found in Data/Modules/ - each text file contains a list of indexes per category/module

2) ConnectomeAnalysis.ipynb:
Jupyter notebook containing relevant code to perform community louvain on connectome-level data, devide connections to categories (subcortical-cortical, interhemispheric, intrahemispheric and intramodular) and calculate connection strength and streamline lenght per connection type

3) Functions.py:
Contains useful functions for connectome analysis.

4) ConnectomeQC&Preprocessing.ipynb:
Notebook to quickly visualise and QC connectomes. 

Feel free to use any of the data and/or code; please consider citing our paper if you do so.
If you have any questions regarding this repo not covered in our paper, please email a.zarkali@ucl.ac.uk. 
