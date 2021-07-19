# FasterHETrace

This repository provides two sets of codes: 

  1) Implementation of faster trace (C++)

	   1-1) bench_unroll.cpp  implementation of our unrolled trace

	   1-2) bench_seq.cpp  implementation of the sequential trace (rotations-and-sums method)

  2) Cost analysis in terms of # modular multiplications (python3)

# 1) Implementation with PALISADE (on the CKKS scheme)

## Prerequisities

  - [Modified version of PALSIADE v1.9.2](https://github.com/uishi/Modified_PALISADEv1.9.2) 

    -> added faster automorphism via simple permutation (`Permute` function in poly.h/cpp and DCRTPoly.h). 
		
		Note that the faster automorphism is applied to both 1-1) and 1-2) to conduct fair comparison.

## Build and run

```
   mkdir build && cd build && cmake .. & make
```

#### Unrolling

      ./bench_unroll [logN] [M] [h] [L] [ell] [Delta] [d] [T]

e.g.

      ./bench_unroll 14 14 7 3 2 35 4 11

#### Rotations-and-sums

     ./bench_seq     [logN] [M]     [L] [ell] [Delta] [d] [T]

e.g.

     ./bench_seq 14 14 3 2 35 4 11

#### Description of each argument

	 log N:  ring dimension

	 M: # iterations for rotations-and-sums

	 h: # unrolled-iterations

	 L: maximal level (the number of RNS moduli for this level is L + 1)

	 ell: level to examine (the number of RNS moduli is \ell + 1 where \ell <= L)

	 Delta: scaling factor for the CKKS encoding

	 d: maximal number of digit for key-switching (an integer that divides L)
	 
	 ---> # of special moduli is determined by k = (L+1)/d

	 T: the number of experiments to perform 
	

# 2) Cost Analysis

  The cost analysis is written in python3.


## Reproducing Figures in Theoretical Analysis

      ./init_resulting_directory.sh

      cd ./cost_analysis
  

### Figure 1 (Rate)

     python3 run_plot_fig1_ell_vs_rate.py

 -> choose 1 and press enter

 Plots are stored in `/tmp/result/no_last_digit_rate/`.


### Figures 2 and 3 (Cost and Its Breakdown)

     python3 run_plot_fig2_fig3_dominant.py

 Plots are stored in `/tmp/result/cost_plot/` and `/tmp/result/breakdown/`.

### Figure 4 (# Eval Keys)

     python3 run_plot_fig4_storage_various_M.py

 Plots are stored in `/tmp/result/numkey/`.
