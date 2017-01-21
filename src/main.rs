#![allow(non_snake_case)]

extern crate rand;
use rand::Rng;
use rand::distributions::{IndependentSample, Range};
mod energy;
use energy::*;
use std::io::prelude::*;

const LJ_EPS : f64 = 1.0;
const LJ_SIG : f64 = 1.0;

const TAILCORR : bool = true;
const SHIFT: bool = false;

// easy printing to stderr
macro_rules! println_stderr(
    ($($arg:tt)*) => { {
        let r = writeln!(&mut ::std::io::stderr(), $($arg)*);
        r.expect("failed printing to stderr");
    } }
);

fn main() {

    // define all the stuff
    let minim_steps  = 1000000;
    let sample_steps = 100000;

    let num_particles: usize = 512;
    let density = 0.7;
    let temperature = 0.9;

    let cutoff = 3.0;
    let displacement = 0.1;


    println_stderr!("");
    println_stderr!("################################################################");
    println_stderr!("##################  LJ Monte Carlo Simulation  #################");
    println_stderr!("################################################################");
    println_stderr!("");


    // initialize stuff
    let beta = 1.0/temperature;
    let volume = (num_particles as f64)/ density;
    let length  = volume.cbrt();
    let (l_x, l_y, l_z) = (length, length, length);
    let cutoff_squared = cutoff * cutoff;

    let mut rng = rand::thread_rng();
    let particle_range = Range::new(0, num_particles-1);

    let mut rx : Vec<f64> = vec![];
    let mut ry : Vec<f64> = vec![];
    let mut rz : Vec<f64> = vec![];
    loop {
        rx.push(l_x * rng.gen::<f64>());
        ry.push(l_y * rng.gen::<f64>());
        rz.push(l_z * rng.gen::<f64>());
        if rx.len() == num_particles { break; }
    }

    let e_corr = if TAILCORR { 8.0/3.0*std::f64::consts::PI*density*LJ_EPS*LJ_SIG.powi(3)*((1.0/3.0*(LJ_SIG/cutoff).powi(9)) - (LJ_SIG/cutoff).powi(3)) } else { 0.0 };
    let p_corr = if TAILCORR { 16.0/3.0*std::f64::consts::PI*density.powi(2)*LJ_EPS*LJ_SIG.powi(3)*((2.0/3.0*(LJ_SIG/cutoff).powi(9)) - (LJ_SIG/cutoff).powi(3)) } else { 0.0 };

    println_stderr!("Particles: {}, Density: {}, Temperature: {}", num_particles, density, temperature);
    println_stderr!("System volume: {:8.3}, Dimensions {:.3}/{:.3}/{:.3}", volume, l_x, l_y, l_z);
    println_stderr!("Minimization steps: {}, Sampling steps: {}", minim_steps, sample_steps);
    println_stderr!("LJ params eps: {}, sigma: {}, cutoff: {}", LJ_EPS, LJ_SIG, cutoff);
    println_stderr!("Tailcorr: {:8.3}, Shift: {:8.3}, Pressurecprr: {:8.3}", e_corr, SHIFT, p_corr);

    let (mut energy, mut virial) = get_total_energy(&rx, &ry, &rz, num_particles, l_x, l_y, l_z, cutoff_squared, e_corr);
    let mut energy_sum = 0.0;
    let mut virial_sum = 0.0;
    let mut step_counter = 0;
    let mut accept_counter = 0;

    println_stderr!("");
    println_stderr!("################################################################");
    println_stderr!("#####################  Energy Minimization  ####################");
    println_stderr!("################################################################");
    println_stderr!("");

    for step in 0..minim_steps+sample_steps {

        // select rnd particle
        let rnd_index = particle_range.ind_sample(&mut rng);

        // store old position
        let oldX = rx[rnd_index];
        let oldY = ry[rnd_index];
        let oldZ = rz[rnd_index];

        // old particle energy
        let (old_particle_energy, old_particle_virial) = get_particle_energy(&rx, &ry, &rz, rnd_index, num_particles, l_x, l_y, l_z, cutoff_squared);

        // rnd displacement and PBC
        rx[rnd_index] += ( rng.gen::<f64>() - 0.5 ) * displacement;
        ry[rnd_index] += ( rng.gen::<f64>() - 0.5 ) * displacement;
        rz[rnd_index] += ( rng.gen::<f64>() - 0.5 ) * displacement;
        if rx[rnd_index] < 0.0 { rx[rnd_index] += l_x }
        if rx[rnd_index] > l_x { rx[rnd_index] -= l_x }
        if ry[rnd_index] < 0.0 { ry[rnd_index] += l_y }
        if ry[rnd_index] > l_y { ry[rnd_index] -= l_y }
        if rz[rnd_index] < 0.0 { rz[rnd_index] += l_z }
        if rz[rnd_index] > l_z { rz[rnd_index] -= l_z }

        // calculate energy difference
        let (new_particle_energy, new_particle_virial) = get_particle_energy(&rx, &ry, &rz, rnd_index, num_particles, l_x, l_y, l_z, cutoff_squared);
        let dE = new_particle_energy - old_particle_energy;

        //accept move
        if rng.gen::<f64>() < (-beta * dE).exp() {
            accept_counter += 1;
            if step % 1000 == 0 { // calculate total energy every 1000 steps to account for rounding errors
                let (e, v) = get_total_energy(&rx, &ry, &rz, num_particles, l_x, l_y, l_z, cutoff_squared, e_corr);
                energy = e;
                virial = v;
            } else {
                energy += dE;
                virial += new_particle_virial - old_particle_virial;
            }
        } else { // or restore old position
            rx[rnd_index] = oldX;
            ry[rnd_index] = oldY;
            rz[rnd_index] = oldZ;
        }

        // update sums for averaging
        step_counter += 1;
        energy_sum += energy;
        virial_sum += virial;


        if step_counter % 5000 == 0 {
            println!("Minim {}\tEnergy: {:.3}\tVirial: {:.3}\tAcceptance:{:.1}\tDisplacement: {:.3}", step_counter, energy, virial, 666, displacement);
        }

        // reset sums for sampling
        if step == minim_steps-1 {
            println!("Starting averaging!");
            step_counter = 0;
            energy_sum = 0.0;
            virial_sum = 0.0;
        }

    }


    let final_energy = energy_sum/step_counter as f64;
    let particle_energy = final_energy / num_particles as f64;
    let final_virial = virial_sum / 3.0 / step_counter as f64 / volume;
    let pressure = virial_sum / 3.0 / step_counter as f64 / volume + density * temperature + p_corr;
    println!("Steps: {}", step_counter );
    println!("Avg Energy: {:.3}", final_energy);
    println!("Energy/Particle: {:.3}", particle_energy);
    println!("Virial: {:.3}", final_virial);
    println!("Pressure: {:.3}", pressure);


}



