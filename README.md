# list-update

This repository contains the code used for the master's thesis "Analysis and Implementation of Algorithms for Self-organizing Lists in the Paid-exchange Model"

The project list-update-final contains the project with the code for the experimental costs of the online algorithms "COUNTER deterministic", "COUNTER randomized", RANDOM RESET and TIMESTAMP.
It also contains the code for the experimental cost difference between the optimum offline algorithm and the pairwise approximation algorithm.

Inside the list-update-final folder there needs to be a folder named ```dataset``` with the corpora to be used in the experiments.

The setup to be used in the experiments can be changed in the ```main.cpp``` file. It can be chosen among static list configuration and dynamic list configuration.
Optionally, it can be choosen to transform the files of the corpora using Burrows-Wheeler transformation.
Similarly the experiments to be run can be changed accordingly. It can be chosen among the experiments calculating the experimental costs and the experiments calculating the cost difference between the optimum offline algorithm and the pairwise approximation algorithm.

This repository also contains three utility jupyter notebooks used to analyze and plot the results collected from the experiments.

If you need help running the experiments and/or the scripts please contact me at anisa dot llaveshi at tum dot de
