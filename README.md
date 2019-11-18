This is the C code used in "The Free and Freer XY model" by Petter Holme and Yerali Gandica. It is an Excange Monte Carlo algorithm using the Metropolis updating. To compile it, first create a directory o for the object files, create directories for the output, and then compile it. 

`mkdir o` 
`mkdir out` for the free XY model 
`mkdir outer` for the freer XY model 
`mkdir out{,er}/{0,1,2,3,4,5,6,7,8,9}{0,1,2,3,4,5,6,7,8,9}` one for every temperature step 
`mkdir out/{0,1,2,3,4,5,6,7,8,9}{0,1,2,3,4,5,6,7,8,9}/{1.0,2.0,4.0,8.0}` one for every avg. degree 
`make` 

The free XY model program takes three arguments: the number of spins, the average degree, a random seed. You can run the freer XY model program in the same way, but without the average degree. For example

`./freexy 64 8.0 2775582768118365280`

The seed should be a 64-bit unsigned integer, and preferably generated by a RNG (i.e. it should be as entropic as possible, and different for every run in simulations).

Since the Metropolis condition mostly serves as a guarantee for the statistics of the sampling, Monte Carlo programs are to some extent more of an art than a science. One is allowed a number of tricks to speed up the code. One such somewhat arbitary trick is the Fermi function to get trial angles. This is done to get the acceptance ratio away from 0 and 1 (when the update dynamics are fastest). Some other tricks and techniques worth noticing:
- For very hight temperatures it is important to design updates such that they sample maximally random configuration. A classic mistake is to have one trial step that undos the previous one (try adding link (i,j), then try deleting (i,j)).
- We represent the angles as 16-bit integers. This allows fast arithmetics (since it is modulo 16-bit, just like the angles are modulo 2 pi), and one can tabulize trigonometric values of 2^16 values.
- We use 1 - cos for the energy rather than -cos. This is in order to not lose precision for very small angles (which are the most important).
- For calculating quantities, we use long double precision. This is because the susceptibility and specific heat requires the subtraction of similar numbers which involves the loss of precision.
- The advantage with exchange Monte Carlo over other parallell tempering schemes is that it guarantees that the Boltzmann distribution is correctly sampled after every temperature swap. The disadvantage is that gaps in the energy could appear that blocks the diffusion of the replicas in the temperature space. We do two things to counteract this: First, we thermalize 10^4 sweeps before the first exchange trial. Second, we have the criteria that we consider averages independent only when they traversed at least a quarter of the temperature levels.
- Our temperatures levels are set in freexy.h. They follow an exponential series determined by the parameters `NT` (number of temperature steps), `TMIN` and `TMAX`.
