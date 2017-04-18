# 2017-02-19 13:14:14
- Local field potential can be inferred from LFPs; [Rasch et al., 2008]
- Origin of local field potential; [Buzsaki et al., 2012]
- spike-LFP coherence; [Mechanisms of Gamma Oscillations]
- {Tried calculate power spectra of LFP}: it seems that all LFP signal in my network stays at a noise level
***
# 2017-02-21 13:14:10
- Rename all existing files; [done]
- Construct a path dependent network program; run and check the result; [not yet]
- write the current situation in my project proposal; [not yet]
- summarize the results got since last time met with Douglas, make a ppt; [not yet]

***
# 2017-02-22 15:27:33
Review about interview with professor Jeffrey Erlich:
	What kind of dynamical pattern that my network can generate;
Review about interview with professor XJWang: [work done by XiaoJing's group]
	1. network dynamic analysis;
	2. cross level neural mechanism understanding;
	3. network modeling in a whole brain region;

***
# 2017-02-24 21:40:00
Talked with Douglas in detail:
1. Differences between passive research and active research; Try to show the strength of my own to defend other's question;
2. As for the experimental setup in LFP, the recording was operating on rats, mainly in hippocampus (CA1); As for previous cases for data analysis, the results are not quite promissing for presenting the relationship between neurons and properties of network; In this case, our group came up with ideas to adapt time-delayed mutual information to analyze LFP data. My work is to test whether this proposal is reasonable and reliable, and eventually, give some distinct ideas for understanding TDMI signal;
3. As my previous work, both network simulating program and TDMI calculation program are worked;[those are done there]
4. Go to check several cases that seem strange previously; including:
	4.1. Estimate the time difference between the time of the peak in TDMI plot and zero time-delayed point. Find its dependence from the parameters of neurons and network;
	4.2. Find the qualitative rules for the relative value of MI with different delayed time-step, as well as the number of histogram's bin;
	4.3. Run simulation based on larger sized network. Analyze the different contribution to LFP from excitatory and inhibitory neurons;
	4.4. Recalculate the TDMI in a specific setup; [in the case of #65 neuron, second-order connection, 1/32ms timing-step;]
	5. Focus on the simplest cases, based on current condition, find the mechanism of MI signal as result of membrane current;
Finished thesis proposal;
***
# 2017-02-25 00:40:37 [Plan]
Upload the thesis proposal; [done]
Send email to Douglas and David; [done]
Answer questions above;
Modify mi.h, mi.cpp and mi_test.cpp
***
# 2017-02-26
Finished mi.h mi.cpp's modification;
***
# 2017-02-27 10:53:48
Start to write summary of previous work, as well as answer the previous questions;[unfinished]
replant all of my program into liunx;
***
# 2017-02-28 10:24:26
configure the plotting program;
plot all the figure used in next report;
plot mi_test.cpp's result;
	1. TDMI of two random Poisson seqence;
	2. TDMI of two double Gaussion random variables; [equally binned, bin size = 0.1]
	3. TDMI of two double Gaussion random variables; [Uniformly occupyed, expected occupance = 10]
	4. TDMI of spikes and LFP from # 65 neruon and its directly connected post-synaptic neurons, bin size 0.01]
	5. TDMI of spikes and LFP form # 65 neuron and its directly connected post-synaptic neurons, expected occupancy 10]
***
# 2017-03-01 10:44:34
rewrite LFP program and program-data interface in C++;
***
# 2017-03-03 00:44:08
Finished lfp in c++;
the efficiency of this program need to be tested;
Upload cods to github;
***
# 2017-03-03 11:13:42
accomplished debugging: lfp.cpp main_lfp.cpp mian_mi.cpp;
***
# 2017-03-03 17:07:26
checking the correctness of new program;
old version finished;
***
# 2017-03-03 20:57:39
It appears that they cannot generate the same answer; 
it migth lead from different random seeds;
Go to Matlab to check the histogram of lfp;
***
# 2017-03-03 16:11:58
There is a huge difference between histogram between these two cases; 
My suspection is that bugs appears at the process of the generation of Poisson feedforward signals;
Tomorrow, first go and check this part;
***
# 2017-03-04 09:16:09
Current test shows that there is no problem in external Poisson Generator;
***
# 2017-03-04 12:14:14
No problem in single neuron program;
Something wrong at neural network program;
They both obey fourth order convergence, while they respond differently for same input;

***

# 2017-03-04 20:13:59
* Current results:
	Comparing those two codes:
	1. as for numerical convergence, all cases obey 4th order numerical convergence, including single neuron case and neural networks, for both old and new codes;
	2. for single neuron case, results between old and new codes are identical;
	3. for neural networks, the results diverges at some cases;
* Current result:
	1. 100 neuron simulation with identical input and network structure: pure excitatory driven; 80% neuron are excitatory; 10% connecting probability between two loop;
	2. simulation time period: 2000 ms;
	3. timing step 2^-5 ms;
	4. #65 neurons connected with 10 neurons: 6 are excitatory, 4 are inhibitory;

***

# 2017-03-05 10:29:55
**What i am going to do this morning:**
* Apply new codes first;
1. 100 neuron; 80% excitatory; 10% interloop connectivity; no rewiring; (neuron type seed 1, 2; Poisson seed, 5, 6, interloop seed 100);
2. 20000 ms simulation time; 2^-5 ms timing step;
3. investigate #65 neuron;
	I. first order connected LFP; (to investigate its physical time delay for the peak of MI)[done]
	II. pure excitatory neuron in 1st order connected cluster;[done]
	III. pure inhibitory;[done]
	IV. 2nd order connection;[done]
	V. record rules for MI magnitude with different timing step and expected occupancy;[done]
	VI. analysis the nonlinear effects that excitory and inhibitory neuron attribute to LFP; [MATLAB]
*NOTE: test old code with stronger input, test new code with less input;*
*Data has been generated, their raster plots are ready to be checked;*

***

# 2017-03-09 14:49:12
* Review about what discussion with Douglas;
	> Test and compare old code and new code again to find the fact which leads to their inconsistancy; [DONE]
	> Test and verify the physical meaning of the time delay between spike and local field potential;	
	> Test the behavior of inhibitory neuronal group: including different pairs of neurons, different network structures;
	> Test different dynamical region: different pattern of inputs as well as different portion between feedforward and interneuron correlation;
* Write configuration file for two network system simulation;[Done]

***

# 2017-03-10 10:58:00
	Test and verify the physical meaning of the time delay between spike and local field potential;	

**Results**

1. different expected occupancy would impact the magnitude of mutual information.
	a. for smaller occupancy, the magnitude of MI becomes larger;
	b. Within certain interval, results based on larger occupancy would have larger signal-to-noise ratio;
	c. In my future test, i would choose 50 as the standard expected occupancy;

2. For current configurations: 

####[PreNet Section]  
	PreNetNeuronNumber = 100  
	PreNetConnectingDensity = 3  
	PreNetTypeProbability = 0.8
	PreNetTypeSeed = 1
	PreNetDrivingType = true
	PreNetDrivingRate = 1.5
	PreNetExternalDrivingSeed = 3
####[PostNet Section]
	PostNetNeuronNumber = 100
	PostNetConnectingDensity = 3
	PostNetTypeProbability = 0.8
	PostNetTypeSeed = 2
	PostNetDrivingType = true
	PostNetDrivingRate = 1.5
	PostNetExternalDrivingSeed = 4
####[Internetwork Section]
	ConnectingProbability = 0.1
	ConnectingSeed = 100
####[Time Section]
	MaximumTime = 10000
	TimingStep = 0.03125

The time delay of the peak of MI is approximately 0.07-0.08 ms. Basically, it was a invariant value under the change of timing-step of TDMI.

3. Based on previous configure, test LFP generated by excitatory neurons;
> for LFP generated by excitatory neurons, the time delay is approximately 0.03ms. Meanwhile, the MI would maintain it peak value for a short period for 0.16 ms;  

***

# 2017-03-11 10:41:34
	Test and verify the physical meaning of the time delay between spike and local field potential;	

**Results**

4. Based on previous configure, the magnitude of MI of spikes and LFP generated by inhibitory neurons are similar to those of spikes and LFP generated by excitatory neuron;
5. For LFP generated by second order connected neurons(71 neurons) of #65 neurons, the MI signal lies on the noise level; for its excitatory portion (56 neurons), things go similarly; for its inhibitory neurons(15 neurons), there is a non significant peak. 
	- The primary guess and conclusion is that excitatory and inhibitory neuron contribute the amount of effect of MI.
	- The magnitude of MI peak is positively related to the number of neurons in current LFP model.[Linear combinition or not, which need derivation.]
	- Second order connected neuron might have much less effect on MI peak.[To be verified]
*To do: Add neuron selection function to lfp.h and lfp.cpp.*
6. For inbibitory neuron, in loop 1, TDMI cannot indicate the relationship between spikes and LFP;

## Plan for tomorrow:
Build Python script to verify those three foundings;
Build neuronal selection functions;
Start to find RA position;

***

# 2017-03-14 19:59:19

- Basically finished python scripts for automatically execute simulation program and output figures;

***

# 2017-03-15 20:42:14

- Fixed bugs of directory in `./multi-network/main.cpp` `./lfp/main.cpp` `./tdmi/main.cpp`
- Add complete version of python scripts for verifying previous three proposals

***

# 2017-03-19 15:58:13

Review about discussion with Douglas:

1. Is the difference between TDMI graph a result from the different dynamical state of pre-network neurons (To see whether there is any different firing state)?
2. The different behavior of post-synaptic neurons that directly connected to the same neuron;
3. different component, different dynamic processes;
4. To analyze the correlation between second order connected neurons;

# 2017-03-21 01:06:04

Finished Report of Week 3 & 4

# 2017-03-21 22:10:40
## Finished
- Fix the old syntax of STL vector with new ones, which increase the efficiency of simulation of two-network system;
- Update the method to calculate LFP by directly input the total membrane current at the begining;

## To do:
1. Define parameters of neuronal dynamics:
	- Firing frequency;
	- connecting neighborhood;
	- Driving component;
2. Define a pattern for TDMI data;
	- noise level [mean value and std];
	- maximum peak value; [difference with noise level]
	- time to reach the peak;
	- time constant for decay side;

# 2017-03-28 11:47:37

- Built program for signal-noise ratio analysis;

## To Do:
1. Run test based on following network structures:
	- two network are both regular;
	- post network are rewired with 0.1 rewiring probability;
	- pre network are rewired with 0.1 rewiring probability;
	- two network are both rewired with 0.1 probability;
2. Labeled each neuron with the clustering coefficienct as well as mean path of its belonged sub neuron clusters. Investigate the change of signal-noise ratio respect to different network structures;
3. Start to write the introductory part of my thesis paper;

# 2017-04-05 10:53:04

## Review for previous progress

1. Run programs for 200 neurons in three different cases:
    - Both networks are regular;
    - Both networks are rewired;
    - Pre network is rewired;
    - Post network is rewired;
2. Data analysis;
    - calculate the mutual information for one degree LFP;
    - record the signal-noise ratio for each case; (signal represents for the maximum value of MI peak)
    - calculate the clustering coefficient and the mean path of the neuorns in their own network as well as across networks; **The basic found is that in the current simulation test, neurons contribute to the local field potential have no self clustered effect, which means that each neuron locate in its own neuron cluster in post network**
    - One strange problem is that neurons with indices close to 200 all have low signal-noise ratio;
3. Orignal thoughts:
    1. Comparing neural signal with identical initial condition and different network structure, to see whether the effect of MI are dominant by cross-network signals or localized dynamics;
    2. To find parameters that can quantify the network structures;
4. TO DO:
    1. why two neurons in the post network connected with the same neuron in pre network have completely different behavior in TDMI signal; *To analysis the phase different between the potential and firing sequency of these two neurons*
    2. Not all excitatory neurons have effective MI peak with their one-degree LFP,

## Done:
	# Standard simulation configuration;
	# [PreNet Section]
	PreNetNeuronNumber = 200
	PreNetConnectingDensity = 6
	PreNetTypeProbability = 0.8
	PreNetTypeSeed = 1
	PreNetDrivingType = true
	PreNetDrivingRate = 1.5
	PreNetExternalDrivingSeed = 3
	PreNetRewiringProbability = 0
	PreNetRewiringSeed = 5
	# [PostNet Section]
	PostNetNeuronNumber = 200
	PostNetConnectingDensity = 6
	PostNetTypeProbability = 0.8
	PostNetTypeSeed = 2
	PostNetDrivingType = true
	PostNetDrivingRate = 1.5
	PostNetExternalDrivingSeed = 4
	PostNetRewiringProbability = 0
	PostNetRewiringSeed = 6
	# [Internetwork Section]
	ConnectingProbability = 0.1
	ConnectingSeed = 100
	# [Time Section]
	MaximumTime = 10000
	TimingStep = 0.03125
1. Run simulation with Poisson driving seed is 30 (for pre network) and 40 (for post network);
2. Run simulation with cross-network configuration genrerating seed 200;

*The results shows that the pattern of maximum signal-noise is independent to initial condition. On the other hand, the pattern of signal-noise ratio for the whole network has a spectial pattern.*

# 2017-04-07 14:26:07

**More Tests**

1. simulation based on network with pure excitatory neurons;
2. Fixed the bug that the program cannot adapte the right order of synaptic input, and run test based on all-excitatory-neuron networks;

# 2017-04-18 10:59:22

1. Run new test, which is stored in .../ResearchData/Apr18/test1/
2. Write python program to anaylze new data;
