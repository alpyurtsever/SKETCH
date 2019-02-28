First run
	CreateData_PhaseRetrieval.m
	CreateData_MaxCut.m
	CreateData_Weather.m
	CreateData_StreamVel.m
	CreateData_NoisyLowRank.m
These scripts draw random instances from the given matrix classes and save
these instances as datasets. Test files will automatically load and use them once the datasets are created. 

There are two alternative ways to run the tests. To run them in cluster using the MATLAB parallel computing toolbox, 
you can use MAIN_CLUSTER_SYN, MAIN_CLUSTER_DATA and MAIN_CLUSTER_ERREST scripts. For this, you should copy the complete 
set of codes, including the datasets into a folder in your cluster. You should also modify the line 43 & 56 & 35 of 
MAIN_CLUSTER_DATA & MAIN_CLUSTER_SYNT & MAIN_CLUSTER_ERREST accordingly. Once you run this script, it will batch the experiments 
using MATLAB SLURM interface. Once the experiments are finished, download the results to your local machine.

If you do not have access to a cluster, you can run the tests in your local 
machine, however this can take some time. You might want to decrease the number of Monte Carlo simulations for this case. 

Once the experiments are finished, you can generate the figures in the paper 
by running files "Plot_Figure_*.m".

To get the results in Fig.3, Fig.7, Fig.8, and FigSM20-25, you need to run
Test_Plot_Figure_7_8
These tests directly run in the local machine, and plots the results. 

Tests for Fig.3 and Fig.SM20 to Fig.SM25 are omitted due to Copyrights of the dataset. 

This toolbox uses some scripts copied (and possibly modified) from PRACTICALSKETCHING toolbox developed in [TYUC2017P].
This toolbox uses some scripts copied (and possibly modified) from Nys-SKETCH toolbox developed in [TYUC2017Nys].
	
NOTE: You should run and save results for the tests only once. When there are 
more than one set of saved results in results folder, plot scripts might 
return an error. In this case, delete the old results for each set of experiments 
and rerun the plot scripts.

See our reference papers for more details.

[TYUC2017P] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher, Practical Sketching Algorithms for Low-Rank Matrix Approximation, Accepted to SIAM J. Matrix Anal. Appl., August 2017.

[TYUC2017Nys] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher. Fixed-Rank Approximation of a Positive-Semidefinite Matrix from Streaming Data. In Proc. 31st Conference on Neural Information Processing Systems (NIPS), Long Beach, CA, USA, December 2017.

[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher. Streaming Low-Rank Matrix Approximation with an Application to Scientific Simulation.

Coded by: Alp Yurtsever
Ecole Polytechnique Federale de Lausanne, Switzerland.
Laboratory for Information and Inference Systems, LIONS.
contact: alp.yurtsever@epfl.ch
Created: August 29, 2016
Last modified: February 22, 2019