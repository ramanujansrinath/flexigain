# Code and data repository for Srinath et al., 2024

## Coordinated Response Modulations Enable Flexible Use of Visual Information

### Abstract
We use sensory information in remarkably flexible ways. We can generalize by ignoring task-irrelevant features, report different features of a stimulus, and use different actions to report a perceptual judgment. These forms of flexible behavior are associated with small modulations of the responses of sensory neurons. While the existence of these response modulations is indisputable, efforts to understand their function have been largely relegated to theory, where they have been posited to change information coding or enable downstream neurons to read out different visual and cognitive information using flexible weights. Here, we tested these ideas using a rich, flexible behavioral paradigm, multi-neuron, multi-area recordings in primary visual cortex (V1) and mid-level visual area V4. We discovered that those response modulations in V4 (but not V1) contain the ingredients necessary to enable flexible behavior, but not via those previously hypothesized mechanisms. Instead, we demonstrated that these response modulations are precisely coordinated across the population such that downstream neurons have ready access to the correct information to flexibly guide behavior without making changes to information coding or synapses. Our results suggest a novel computational role for task-dependent response modulations: they enable flexible behavior by changing the information that gets out of a sensory area, not by changing information coding within it. 

**Note:** Most of these functions don't just generate the figures with pre-analyzed data; wherever possible, I have uploaded the raw data and analysis code. Please [email me](mailto:ramsrinath@uchicago.edu) if you find missing dependencies or need access to even raw-er data like spike times or LFP data. For some analyses, the data was too large to fit in this repository. I will upload those on OSF (and the rest of the `data` folder) and link here.

### Code Components
- [x] Continuous curvature estimation task
    - [x] Shape generation (requires [geom2d toolbox](https://www.mathworks.com/matlabcentral/fileexchange/7844-geom2d))
    - [x] Behavioral analysis
    - [x] Curvature tuning and selectivity
    - [x] PCA visualization and curvature decoding
    - [x] Out-of-set or across-shape curvature decoding
    - [x] Comparison of shape-specific and shape-general decoding
    - [x] Saccade decoding

- [x] RNN modeling
    - [x] RNN training python code (requires [PsychRNN](https://psychrnn.readthedocs.io/en/latest/) and dependencies)
        - [x] Simple inputs
        - [x] CNN inputs
    - [x] Decoding curvature with VGG16 activations
        - [x] Random shape training
        - [x] Exhaustive shape training
    - [x] Analysis
        - [x] Arc-dependent curvature report for single shape
        - [x] Arc-dependent curvature report with VGG-16 inputs
            - [x] Random shape training
            - [x] Exhaustive shape training
            - Note: The data for this analysis will be posted elsewhere, probably OSF. I'll update this and add a link when it's ready. Meanwhile, you can email me and I can share it.

- [x] Curvature-Color task
    - [x] Stimulus generation
    - [x] Psychophysics
    - [x] Neural responses to each shape
    - [x] Gain vs selectivity
    - [x] PCA/QR visualization
    - [x] Curvature/color decoding and relationship to choice

- [x] Tuning and fixed readout simulations
    - [x] Arc-dependent curvature report across shapes
    - [x] Feature-dependent curvature/color 2AFC task
