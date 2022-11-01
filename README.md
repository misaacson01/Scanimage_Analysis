# Scanimage_Analysis
 
After imaging has been collected using the SI_automated script, run these data processing steps in order:

1. PreProcess_Part1 (motion correction, spectral unmixing)
2. Suite2p (neuron segmentation, dF/F and spike extraction)
3. (outdated) PreProcess_Part2 (Creating an image stack to help link neurons between timepoints and IDing cell types)
4. cellLinker_Part3 (GUI to link neurons between timepoints, ID cell types)
5. (TBD) Analysis_Part4: (plot spontaneous activity levels, orientation/direction tuning, individual cell changes across timepoints in a single mouse)
6. (TBD) Analysis_Part5: (plot analyzed individual cell data across multiple mice)
