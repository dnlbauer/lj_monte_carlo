# lj_monte_carlo
Simple Lennard Jones Liquid MC simulation implemented in Rust

To compile, run:
```
/ release builds are much faster because of enabled lto
cargo build --release /
```

The compiled binary can then be found in ```target/release.```

# Phase diagram of a 6-12 lj fluid
(cutoff rc of 3.5 sigma)  
![phase diagram of a lj fluid](https://www.researchgate.net/profile/Billy_Todd/publication/7525791/figure/fig1/AS:280682271133696@1443931276144/FIG-1-Phase-diagram-for-the-6-12-Lennard-Jones-fluid-with-a-cutoff-radius-of-r-c-35.png)  
Source: Ge J. et ak., Scaling behavior for a pressure and energy of shearing fluids. Phys Rev 061201, p. 67-68, 2003


# Results
Some results obtained during the creation of this project


![Function of state](results/FoS/FoS.png)
Function of state at different temperatures. 512 particles were simulated with applied potential shift and tailcorrection. Displacement scaled to have 33% acceptance rate during 1mio steps. Averages over 100k steps.

![Minimization](results/energy_minimization/energy_minimization.png)
