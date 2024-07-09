# Code and data repository for Srinath et al., 2024

## Coordinated Response Modulations Enable Flexible Use of Visual Information

### Abstract
We use sensory information in remarkably flexible ways. We can generalize by ignoring task-irrelevant features, report different features of a stimulus, and use different actions to report a perceptual judgment. These forms of flexible behavior are associated with small modulations of the responses of sensory neurons. While the existence of these response modulations is indisputable, efforts to understand their function have been largely relegated to theory, where they have been posited to change information coding or enable downstream neurons to read out different visual and cognitive information using flexible weights. Here, we tested these ideas using a rich, flexible behavioral paradigm, multi-neuron, multi-area recordings in primary visual cortex (V1) and mid-level visual area V4. We discovered that those response modulations in V4 (but not V1) contain the ingredients necessary to enable flexible behavior, but not via those previously hypothesized mechanisms. Instead, we demonstrated that these response modulations are precisely coordinated across the population such that downstream neurons have ready access to the correct information to flexibly guide behavior without making changes to information coding or synapses. Our results suggest a novel computational role for task-dependent response modulations: they enable flexible behavior by changing the information that gets out of a sensory area, not by changing information coding within it. 

### Code Components
- Continuous curvature estimation task
    - Shape generation (requires geom2d toolbox)
    - Behavioral analysis
    - Curvature tuning and selectivity
    - PCA visualization and curvature decoding
    - Out-of-set or across-shape curvature decoding
    - Comparison of shape-specific and shape-general decoding
    - Saccade decoding

- RNN modeling
    - Arc-dependent curvature report for single shape
    - Arc-dependent curvature report with VGG-16 inputs
        - Random shape training
        - Exhaustive shape training

- Curvature-Color task
    - Stimulus generation
    - Psychophysics
    - Feature tuning and selectivity
    - PCA/QR visualization
    - Curvature/color decoding and relationship to choice

- Tuning and fixed readout simulations
    - Arc-dependent curvature report across shapes
    - Feature-dependent curvature/color 2AFC task
