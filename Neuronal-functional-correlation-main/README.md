# Neuronal-functional-correlation

Code and Data for "Stable visual representations correlate with mesoscale neuronal coordination"

Project Overview

This repository provides the source code and analysis pipeline for our study on how the brain maintains stable representations despite neural variability. We address the fundamental question of how stable perception arises from noisy, single-neuron activity.

Our work proposes and validates that the solution lies not in individual neuron activity, but in **dynamic neuronal functional correlation (NFC)**—the transient, coordinated interactions between neurons. This repository contains the tools to replicate our findings, which demonstrate that these network-level patterns provide a robust and stable neural code.


Directory Structure

*   **Example Data**: This folder contains a subset of the data used in the study, provided as examples to illustrate the data format and structure.

*   **Dynamic NFC Computation**: A framework to calculate dynamic Neuronal Functional Connectivity (NFC) from large-scale neural recordings from RUSH 3D

*   **Stability and Decoding Analysis**: Scripts to quantify trial-to-trial variability and show that NFC-based codes are significantly more stable and yield higher decoding accuracy than codes based on single-neuron activity.

*   **Network Dynamics Characterization**: Tools to analyze the spatiotemporal evolution of functional networks, including the propagation of connectivity waves from visual to higher-order cortices and the characterization of their small-world properties.

*   **In-Silico Modeling**: An implementation of the computational model used to corroborate our experimental findings and test the principles of stable representation through dynamic connectivity.

By providing this code, we aim to facilitate further research into network-centric coding principles and their potential applications in developing more robust artificial intelligence.

For a detailed description of the methodology and findings, please refer to our paper:
> ***

Contributing
We welcome contributions to improve the code and analyses. If you find any issues or have suggestions for enhancements, please open an issue or submit a pull request.
