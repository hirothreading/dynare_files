# Notes on BLL (2013, AER)
My replication codes are basically an update/fix to the BLL replication package which was unable to run on modern Dynare. The main issue is that the version of Dynare that the authors used, Dynare v3.065, allowed the use of matrices in the model block. As far as I understand, this feature was removed from Dynare version 4 onwards. The authors also used modfied versions of ```stoch_simul.m``` and ```dynare_estimation.m``` in order to save IRFs for plotting. I have done my best to update their codes and modifications so that all IRF tensors are contructed and saved in the ```IRF_convert.m``` and ```make_IRFs.m``` files. This means that users can run simulations on modern versions of Dynare without having to modify ```stoch_simul.m``` and ```dynare_estimation.m```.  

Please note that the results from my replication differs slightly from the results in the BLL paper for a few reasons:

1. I have set parameter values for the simulated model (not the estimated model) at their prior means. 
2. I have tried my best to correct any typos found in the BLL replication package. I have left comments throughout the Dynare model files and associated MATLAB ```[model]steady_state.m``` files to record the typos I caught.
3. As some parameters are different due to typo corrections, the estimation strategy (most notably the priors used) is quite different from the original BLL paper.  

To run the Bayesian estimation exercise on your own computer, I suggest you take advantage of the ```parallel``` options available in Dynare. Further information can be found on Willi Mutschler [website](https://mutschler.eu/dynare/installation/parallel/), or ask for advice on the [Dynare forums](https://forum.dynare.org). For reference, the estimation took about three and a half hours on my Apple MacBook Air with an M1 processor and 16GB of RAM, and with macOS Sequoia 15.2, MATLAB r2024b, and Dynare 6.2.

Finally, I have probably made some typos or mistakes of my own. If you spot any, please do email me or reach out to me on the Dynare forums!