# Spindle Iterative RevIsion System (SIRIS)
The Spindle Iterative RevIsion System (SIRIS) is an open source project for human-machine coupled sleep spindle detection by iterative revision. It can be run in a MATLAB environment, and complements the paper at https://doi.org/10.1101/454488.

---

## Graphical User Interface
SIRIS is integrated into an intuitive Graphical User Interface (GUI) for researchers to use conveniently. All fundamental actions that can be executed by SIRIS are accessible from this interface, located in "siris/GUI/Choose_Spindles.m." These include the following:

### Visualization of EEG traces
The GUI for SIRIS contains a simple but effective visualization interface for EEG traces. Simply select *Load Raw Data*, which reads the raw EEG trace **(in .edf format)** into the system. *Load Scores* is optional, and can load timestamps of previous sleep-scored EEG segments **(in .txt format)** into the system. *Plot Graphs* plots the EEG trace on the given figure (the upper axes corresponding to the raw trace and the lower axes corresponding to a 8-18Hz bandpass filtered signal, the filter of which is given in the folder "siris/Helper Functions/BP8_18Hz.mat").

### Superimposement of Spindle Labelings
SIRIS contains, in addition to the visualization of raw EEG traces, a function of loading pre-labeled sleep spindle segments. Select *Load Other Spindles* to load previously marked spindle sets in a read-only format. Select *Load Spindles* to load a previously marked spindle set into the system. Spindle datasets should be **in .m format** in which the timestamps for start and stop times (in frames) of sleep spindles are recorded in a matrix variable named "spindle_points." The spindle sets that the STFT algorithm outputs and those that are labeled through the graphical interface are all automatically formatted to be compatible with this system.

### On-The-Fly Spindle Labeling System
Users are able to easily create new spindle label sets or modify existing sets using the graphical interface. Once the data has been loaded into the system, the user can mark spindles by simply clicking, inside the upper axes, the points for the desired start and end times. Green lines will appear to aid in the process of selecting the precise timestamps, and the user can drag these bars to achieve labeling of the desired section for the sleep spindle. To save the spindle, one can click *Record Spindle* at the top-right corner of the GUI, which appears after the user has marked the desired timestamps in the plot. After this, the timestamps are automatically saved into a file holding all spindle labels from this trial or into the spindle set loaded into the system for editing. After saving each spindle, the timestamps will appear in the interface as dashed black lines. To delete previously labeled spindles, simply select a region in which spindles are to be deleted, and choose *Delete Spindle(s)*. All marked spindles in this region will be deleted, and the corresponding file will be updated.

### Analysis of Labelings
The integrated GUI also supports the function of analyzing the similarity and match of different sets of spindle labelings. After loading the sets to be compared through *Load Other Spindles*, users can simply select *Calculate Statistics*, and choose the desired two spindle sets to be compared. Certain statistics of a by-event analysis, including the F1 score, will be calculated, and a plot of the overlap percentages for spindle events will be shown.

### Batch Processing Using STFT Algorithm
The SIRIS GUI is also able to perform a batch process on raw EEG data using an implemented Short-Time Fourier Transform algorithm (STFT). Simply select *Run Algorithm on Data*, and after inputting a range of STFT parameters to test for, the machine will run in a batch process a sleep spindle detection using specified parameters. Allow the algorithm some time to run on the computer (minutes to hours, depending on the step size and range of thresholds specified).

### Batch Analysis ("Selection") Based on Reference Spindle Set
The SIRIS GUI is able to perform automated selection of STFT parameters given initial (or Revised) spindle labels and the set of all labels obtained from varying STFT parameters. After selecting *Performance Analysis* and inputting desired parameters in the STFT batch processing and the False-Negative tolerance rate, the user can find the optimal STFT parameters that is able to complement the current set. This set can then be used for "Revision."

### Revision
Revision can easily be performed in the SIRIS GUI. After loading two (or several) spindle sets that can contribute to the next reference set as *Other Spindles*, the user can click *Review Labels*, and select those sets that are to undergo Revision. This will bring the user to the Revision interface, where the user is blind of the origin set of each event. For each potential event displayed, the user can either *Accept*, *Reject*, or *Update* the spindle timestamps, forming the corresponding entry in the final output set of Revision.

---

## STFT Algorithm
The STFT algorithm used in SIRIS is described in detail in the complementary paper:

> An automatic spindle detection algorithm was implemented using MATLAB (R2018a, Mathworks, Inc.) based on the Short-Time Fourier Transform (STFT) (Gorur D, Halicil U, Aydin H, Ongun G, Ozgen F, Leblebicioglu K (2002) Sleep spindles detecton using short time fourier transform and neural networks. Neural Networks, 2002. IJCNN ‘02.), achieving fast machine labeling of spindles (∼1 minute of running time for 6 hours of EEG recording). The EEG data was first preprocessed by smoothing and noise reduction. Previously sleep-scored non-REM segments of sleep were then transformed by STFT with a 300ms Hamming window and 250ms overlap. The power of the spindle frequency band (8-16 Hz) was calculated relative to the total power of the signal, and a double-threshold was applied on this power ratio for spindle detection. Segments 0.5-3 seconds long that crossed a lower threshold during their entire lengths and crossed an upper threshold at least once in their durations were considered spindles.
