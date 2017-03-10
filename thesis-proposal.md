#PURPOSES:
##Computational research in neuroscience:
As a fresh member in brain science, computational neuroscience applies information theory and computational tools to help understand neural mechanism and brain functions better. Meanwhile, close relationship between computational neuroscience and experimental neuroscience allows people to develop more computational technology and theories to analyze the vast pool of available neurobiology data and together to build a view of how brain works in a whole-brain level.

##Integrate and fire neuron model:
As the simplest version of spiking neuron model, integrate-and-fire neuron model plays a significant role in studies of network dynamics. It eliminates the complex intracellular properties while preserve the important electrophysiological properties of neuron, which tells how neuron process and conduct information. Despite of its limitation to describe the detailed intracellular structure of neuron, it is widely accepted by various work about neural information processing and nonlinear neural network dynamics. As a canonical version of spiking neuronal model, starting from Lapicque's model, integrate-and-fire model diverges into various versions and stimulates wide applications.
In this project, integrate-and-fire neuron model would be applied as the basis of modeling. Mechanism of interneuron communication and computation based on different sparse network structure would be discussed as well as other network dynamics. The details of model in simulation would be modified according to actual experimental needs.

##Local field potential:
In brain data recording, there are several techniques, including electroencephalography(EEG), magnetoencephalography(MEG), electrocorticography(ECoG) and local field potential(LFP). All those methods collect information attributed by extracellular fields and currents and intracellular activities, such as synaptic activities, fast action potentials and intrinsic currents and resonances. Among those, LFP is recorded by small-sized electrodes in the brain. The data combines mainly the low frequency brain rhythm which generate from neurons within approximately 50-350 micrometer from the electrode and slow ionic events from within 0.5-3 mm from the electrode. The potential to achieve high spatial resolution for brain recording becomes a reason to apply this method. On the other hand, compared with superficial layer sampling of EEG, MEG and ECoG, LFP have abilities to record deeper located events in cortex. As a result, its electrophysiological origins are complicated which contains various of components and it would be modeled and discussed in this project in detail.

##Time delayed mutual information:
In information theory, mutual information(MI) is used to quantify the relationship between two sets of data. In principle, the value of MI is invariant under homeomorphic transformation, and comparing with cross-correlation, it has a better performance dealing with data under nonlinear mapping. And in neural information analysis, it is a great tool for identify neuron-neuron correlation based on neuronal signal recording. And somehow, it would be a good tool to investigate the dynamic of single neuron and its related local field potential, which would be covered in this work. In addition, by adding time-delays onto MI analysis, we can observe the change MI in time domain, which usually called time-delay mutual information(TDMI). In this case, we can observe the direction of information flow in different area in neural network. Thus, we can somehow reveal the network structure based on available experimental. And this is the key purpose of this project.

#CONTENT:
1. Understand integrate-and-fire neuronal model and TDMI analysis for neural data;
2. Build programs for IF neuronal network simulation and TDMI analysis;
3. Apply method and program on experimental data;
4. Modified electrophysiological model of local field potential;
5. Verify the effectiveness of TDMI analysis of LFP and show its mechanism;

#METHODS:
1. Theoretical parts:
1.1 Build integrate-and-fire neuronal network simulation programs. Run simulations for sparse connected neuronal network, which is constructed based on the form of 'small-world' network. 
1.1.1 For integrate-and-fire neuronal mode, it is a conductance-based spiking neuronal model with fixed threshold membrane potential and refractory period. Each neuron is driven by both feedforward homogeneous (or nonhomogeneous) Poisson input and interneuron spikes (which could be from either excitatory and inhibitory).
1.1.2 Neuronal network model is arranged in the form of 'small-world' network. Starting from a one dimensional lattice loop, with certain rewiring operation between neurons, networks with certain small-worldness are achieved, which fits the experimental finding about neuronal networks in cortex.
1.2 Build models for local field potential. Use simulating data to test the model;
1.3 Estabish proper program to extimate MI numerically;
2. Experimental parts:
Combining with experimental data set, test the method of TDMI analysis, and compare it with those based on simulating data;



#EXPECTATION:
1. Summarize rules for joint analysis of spikes and local field potential in terms of TDMI;
2. Give detailed explanation to experimental results;
3. Apply TDMI method on experimental data and compare results with those from simulation. Modify the electrophysiological model of local field potential with the help of ;
4. Verify the effectiveness of TDMI analysis of LFP and show its mechanism;
5. Build primary programs that can operate joint analysis of spikes and LFP in terms of TDMI;

Brief summary of purpose of this project:
1. Help understand biophysical origin of local field potential;
2. Give inspirations to experimental LFP data analysis about how TDMI would work within those processes;
3. Find rules for modes correspond to different network structures;
4. Understanding about spiking neuronal model simulation and neuronal data analysis better.

#PLANS:
- Week 1: Accomplish the thesis proposal; Read Papers to understand the electrophysiological mechanism of local field potential in different experimental setups; Define the first version of biophysical model (including its related network structure) and applied it to simulating neuronal data;
- Week 2-8: Test LFP model based on simulating data, and modified it step by step;
- Week 9 (Apr 23): Mid-term thesis defense;
- Week 10-15: Apply TDMI method on experimental data and compare results with those from simulation.
- Week 16 (Jun 7): Final thesis defense;

#Reference：
1. Gerstner, W., & Kistler, W. M. (2002). Spiking neuron models: Single neurons, populations, plasticity. Cambridge university press.
2. Dayan, P., & Abbott, L. F. (2001). Theoretical neuroscience (Vol. 806). Cambridge, MA: MIT Press.
3. Shelley, M. J., & Tao, L. (2001). Efficient and accurate time-stepping schemes for integrate-and-fire neuronal networks. Journal of Computational Neuroscience, 11(2), 111-119.
4. Hansel, D., Mato, G., Meunier, C., & Neltner, L. (1998). On numerical simulations of integrate-and-fire neural networks. Neural Computation, 10(2), 467-483.
5. Rangan, A. V., & Cai, D. (2007). Fast numerical methods for simulating large-scale integrate-and-fire neuronal networks. Journal of Computational Neuroscience, 22(1), 81-100.
6. Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of ‘small-world’networks. nature, 393(6684), 440-442.
7. Cellucci, C. J., Albano, A. M., & Rapp, P. E. (2005). Statistical validation of mutual information calculations: Comparison of alternative numerical algorithms. Physical Review E, 71(6), 066208.
8. Lapicque, L. (1907). Recherches quantitatives sur l’excitation électrique des nerfs traitée comme une polarisation. J. Physiol. Pathol. Gen, 9(1), 620-635.
9. Tuckwell, H. C. (1988). Introduction to theoretical neurobiology. Vol. 1, Linear cable theory and dendritic structure and stochastic theories.
10. Taghva, A., Song, D., Hampson, R. E., Deadwyler, S. A., & Berger, T. W. (2012). Determination of Relevant Neuron–Neuron Connections for Neural Prosthetics Using Time-Delayed Mutual Information: Tutorial and Preliminary Results. World neurosurgery, 78(6), 618-630.
11. Buzsáki, G., Anastassiou, C. A., & Koch, C. (2012). The origin of extracellular fields and currents—EEG, ECoG, LFP and spikes. Nature reviews neuroscience, 13(6), 407-420.
