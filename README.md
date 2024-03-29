# AMP & RAMP

This project implements the AMP and RAMP algorithms for the paper:
- Anonymous Author(s), "Efficient and Effective Algorithms for A Family of Influence
  Maximization Problems with A Matroid Constraint".
  
###This project can be run with one click. Please refer to the following details.

## How to compile

``cmake version >= 3.10`` is required.

``cmake .``; ``make``.

## How to run

Command Format: ``./test_[mrim/rm/aa] [graph_file] [args] [0/1]``

``[mrim/rm/aa]``: The IM-GM instance in which the algorithms will run.

``[graph_file]``: The filename of the graph data (the format can be referred to ``data/facebook.txt``).

``[args]``: The arguments in the IM-GM instance.

``[0/1]``: 
- 0 = Run AMP and competitors (Greedy, Local-Greedy, Threshold-Greedy) by varying #RR sets from $2^1$ to $2^17$ (or from $2^5$ to $2^22$ in RM); 
- 1 = Run RAMP and competitors (RMA / CR-NAIMM / AAIMM) by varying $\epsilon$ from 0.5 to 0.1.

E.g., assume that the graph file is ``data/facebook.txt``:

**Run in MRIM with k=100 and T=20**: ``./test_mrim data/facebook.txt 100 20 [0/1]``

**Run in RM with T=10**: ``./test_rm data/facebook.txt 10 [0/1]``

**Run in ADVIM with k_v=500 and k_e=1000**: ``./test_aa data/facebook.txt 500 1000 [0/1]``

**Run in MRIM with k=2000 (case T=1 in MRIM)**: ``./test_mrim data/facebook.txt 2000 1 [0/1]``

**Note**: In RM, parameters $\alpha_i$ and $k_i$ are fixed to 1.

## Dataset Sources

The datasets used in our paper can be downloaded via the following links. You might need to do extra cleaning for the raw data to align with our input format (see ``data/facebook.txt``).

Facebook: http://snap.stanford.edu/data/ego-Facebook.html

Wiki: http://snap.stanford.edu/data/wiki-Vote.html

Twitter: http://snap.stanford.edu/data/ego-Twitter.html

Google+: http://snap.stanford.edu/data/ego-Gplus.html

Pokec: http://snap.stanford.edu/data/soc-Pokec.html

LiceJournal: http://snap.stanford.edu/data/soc-LiveJournal1.html

We encourage you to try other datasets that can be used in IM tasks to test the performance.