# Faster Homomorphic Trace-type Function Evaluation

This repository provides two sets of codes: 

  1) Implementation of faster trace (C++)

	   1-1) bench_unroll.cpp  implementation of the unrolled trace

	   1-2) bench_seq.cpp  implementation of the sequential trace (rotations-and-sums method)

  2) Cost analysis in terms of # modular multiplications (python3)

# 1) Implementation with PALISADE (using the CKKS scheme)

## Prerequisities

  - [Modified version of PALSIADE v1.9.2](https://fs.yama.info.waseda.ac.jp/yuishi/modified_palisadev1.9.2) 

    -> added faster automorphism via simple permutation (`Permute` function in poly.h/cpp and DCRTPoly.h). 
		
The faster automorphism is applied to both 1-1) and 1-2) to conduct fair comparison.

## Build and run

```
   mkdir build && cd build && cmake .. && make
```

#### Unrolling

      ./bench_unroll [logN] [M] [h] [L] [ell] [Delta] [d] [T]

e.g.

      ./bench_unroll 14 14 7 3 2 35 4 11

#### Rotations-and-sums

     ./bench_seq     [logN] [M]     [L] [ell] [Delta] [d] [T]

e.g.

     ./bench_seq 14 14 3 2 35 4 11

#### Description of each param

	 log N:  ring dimension

	 M: # iterations for rotations-and-sums

	 h: # unrolled-iterations

	 L: maximal level (# of RNS moduli for this level = L + 1)

	 ell: level to examine (# of RNS moduli for this level = \ell + 1 where \ell <= L)

	 Delta: scaling factor for the CKKS encoding

	 d: maximal number of digits for key-switching (an integer that divides L+1)
	 
	 ---> # of special moduli k = (L+1)/d

	 T: the number of experiments to perform 
	

# 2) Cost Analysis 

## Reproducing Figures in Theoretical Analysis

      ./init_resulting_directory.sh

      cd ./cost_analysis
  

### Figure 1 (Rate)

     python3 run_plot_fig1_ell_vs_rate.py

 -> choose 1 and press enter

 Plots are saved in `/tmp/result/no_last_digit_rate/`.


### Figures 2 and 3 (Cost and Its Breakdown)

     python3 run_plot_fig2_fig3_dominant.py

 Plots are saved in `/tmp/result/cost_plot/` and `/tmp/result/breakdown/`.

### Figure 4 (# Eval Keys)

     python3 run_plot_fig4_storage_various_M.py

 Plots are saved in `/tmp/result/numkey/`.
