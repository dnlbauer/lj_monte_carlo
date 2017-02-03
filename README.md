# lj_monte_carlo
Simple Lennard Jones Liquid MC simulation implemented in Rust

To compile, run:
```
// release builds are much faster because of enabled lto
cargo build --release
```

The compiled binary can then be found in ```target/release.```

## Phase diagram of a 6-12 lj fluid
(cutoff rc of 3.5 sigma)  
![phase diagram of a lj fluid](https://www.researchgate.net/profile/Billy_Todd/publication/7525791/figure/fig1/AS:280682271133696@1443931276144/FIG-1-Phase-diagram-for-the-6-12-Lennard-Jones-fluid-with-a-cutoff-radius-of-r-c-35.png)  
Source: Ge J. et ak., Scaling behavior for a pressure and energy of shearing fluids. Phys Rev 061201, p. 67-68, 2003


# Results
Some results obtained during the creation of this project

## Function of State
Function of state at different temperatures. 512 particles were simulated with applied potential shift and tailcorrection. Displacement scaled to have 33% acceptance rate during 1mio steps. Averages over 100k steps.
![Function of state](results/FoS/FoS.png)  

![Minimization](results/energy_minimization/energy_minimization.png)
Energy minimization for a box of 2048 particles at T=0.9 and Density=0.7

## Coexisting liquid/vapor phases
2048 particles were inserted into a box at T=1.0 and density=0.85 (Dimensions 13.406/13.406/40.218), no tailcorrection, cutoff 2.5. The particles filled 1/3 of the box with vacuum space above and below the liquid slab. Finally, a simulation was run for 50000+5000 Cycles (each 2048 steps) writing a frame every cycle.

Density profile avg. over the last 5000 cycles.
![Density profile](results/coexisting/density_z.png)


Chemical potential over all 55000 cycles. z=11.07 and z=27.57 were chosen as liquid boundaries. For every frame, 100 particles insertions per phase were performed.
![Chemical potential](results/coexisting/chem_pot.png)

Surface tension: Running average over the last 5000 Cycles. **TODO shift virial**
![Surface tension](results/coexisting/surface_tension.png)
