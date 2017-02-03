mod trajectory;
use trajectory::*;
mod energy;
use energy::*;
extern crate rand;
use rand::Rng;
use std::env;

const LJ_EPS : f64 = 1.0;
const LJ_SIG : f64 = 1.0;

static MKSA_PLANCKS_CONSTANT_H : f64 = 1.0;
static MASS : f64 = 1.0;

const RUN_AVG_SIZE : usize = 1;

const SHIFT : bool = true;

fn eval_ideal_potential(temperature: f64, lj_eps: f64, volume: f64, particles: f64, thermal_wavelength3: f64) -> f64 {
    let density = volume/particles;
    return -temperature/lj_eps * ( density*thermal_wavelength3 ).ln();
}

#[test]
fn test_eval_ideal_potential() {
    let expected = 12.476649;
    let result = eval_ideal_potential(2.0,1.0,10.0,512.0,0.1);
    assert!((expected-result).abs() < 0.00001, "{}", result);
}

fn main() {
    // parse args
    let args: Vec<String> = env::args().collect();
    let mut filename = "montecarlo.xyz".to_string();
    let mut skip: usize = 0;
    let mut insertions: usize = 100;

    // liquid phase boundaries
    let mut liquid_start = 0.0;
    let mut liquid_end = 0.0;
    for i in 0..args.len() {
        if args[i] == "-f" {
            filename = args[i + 1].clone();
        } else if args[i] == "-s" {
            skip = args[i + 1].parse::<usize>().unwrap();
        } else if args[i] == "-ls" {
            liquid_start = args[i + 1].parse::<f64>().unwrap();
        } else if args[i] == "-le" {
            liquid_end = args[i + 1].parse::<f64>().unwrap();
        } else if args[i] == "-n" {
            insertions = args[i+1].parse::<usize>().unwrap();
        }
    }

    // open file and skip to requiested position
    let mut trj_reader = TrjReader::new(&filename);
    if skip > 0 { trj_reader.skip(skip) };

    // trajectory information
    let mut frame = trj_reader.next_frame();
    println!("{:?}", frame);

    // get some non changing values
    let volume = frame.box_x * frame.box_y * frame.box_z;
    let beta = 1.0/frame.temperature;
    let cutoff_sqr = frame.lj_cutoff * frame.lj_cutoff;

    // rnd number generator for particle insertion
    let mut rng = rand::thread_rng();

    // calculate liquid and gas volume
    let liquid_height : f64 = liquid_end - liquid_start;
    let liquid_volume = liquid_height * frame.box_x * frame.box_y;
    let gas_volume = volume - liquid_volume;

    // thermal wavelength to the power of 3 = h/sqrt(2*pi*m*kb*T) = (2*pi*m*kb*T/h^2)^3/2
    let tw3 =  ((2.0 * std::f64::consts::PI * MASS * frame.temperature / frame.lj_eps)/MKSA_PLANCKS_CONSTANT_H.powi(2)).powf(3.0/2.0);

    // LJ shift
    let e_shift = if SHIFT { 4.0 * LJ_EPS * ( (LJ_SIG/frame.lj_cutoff).powi(12) - (LJ_SIG/frame.lj_cutoff).powi(6) ) } else { 0.0 };


    // average counters
    let mut frame_count = 0;
    let mut widom_sum_gas = 0.0;
    let mut ideal_pot_gas_sum = 0.0;
    let mut widom_sum_liquid = 0.0;
    let mut ideal_pot_liquid_sum = 0.0;
    let mut avg_count = 0;

    // loop over all frames
    let mut counter = 0;
    loop  {
        frame_count += 1;

        // calculate number of particles in each phase
        let mut liquid_count = 0.0;
        let mut gas_count = 0.0;
        for i in 0..frame.num_particles {
            if frame.rz[i] > liquid_start && frame.rz[i] < liquid_end { liquid_count += 1.0; }
            else { gas_count += 1.0; }
        }

        // test particle insertion
        for i in 0..insertions {
            avg_count += 1;

            // liquid test partcile
            let lx = frame.box_x * rng.gen::<f64>();
            let ly = frame.box_y * rng.gen::<f64>();
            let lz = liquid_start + (liquid_height * rng.gen::<f64>());
            let widom_e_liquid = get_particle_insertion_energy(&frame.rx, &frame.ry, &frame.rz, frame.num_particles, lx, ly, lz, frame.box_x, frame.box_y, frame.box_z, cutoff_sqr, e_shift);
            widom_sum_liquid += (-beta*widom_e_liquid).exp();

            // gas test particle
            let gx = frame.box_x * rng.gen::<f64>();
            let gy = frame.box_y * rng.gen::<f64>();
            let mut gz = frame.box_z * rng.gen::<f64>();
            while gz > liquid_start && gz < liquid_end { // retry until we have a particle in gas
                 gz= frame.box_z * rng.gen::<f64>();
            }

            let widom_e_gas = get_particle_insertion_energy(&frame.rx, &frame.ry, &frame.rz, frame.num_particles, gx, gy, gz, frame.box_x, frame.box_y, frame.box_z, cutoff_sqr, e_shift);
            widom_sum_gas += (-beta*widom_e_gas).exp();

            // calculate ideal potentials
            ideal_pot_gas_sum += eval_ideal_potential(frame.temperature, frame.lj_eps, gas_volume, gas_count, tw3);
            ideal_pot_liquid_sum += eval_ideal_potential(frame.temperature, frame.lj_eps, liquid_volume, liquid_count, tw3);

        }

        // print running averages
        if avg_count / insertions > RUN_AVG_SIZE {
            let ideal_gas_potential = ideal_pot_gas_sum / avg_count as f64;
            let ideal_liquid_potential = ideal_pot_liquid_sum / avg_count as f64;
            let excess_gas = -(widom_sum_gas/avg_count as f64).ln()/beta;
            let excess_liquid = -(widom_sum_liquid/avg_count as f64).ln()/beta;
            let mut gas_total = ideal_gas_potential + excess_gas;
            let liquid_total = ideal_liquid_potential + excess_liquid;
            println!("Frame {}\tg_ex: {:5}\tl_ex: {:5}\tg_tot: {:5}\tl_tot: {:5}\t\tnparticles: {}/{}", frame_count, excess_gas, excess_liquid, if gas_total.is_infinite() { 0.0 } else { gas_total } , liquid_total, gas_count, liquid_count);
//            println!("Frame {}\tgas {}\tliquid {}\t particles gas/liquid:{}/{}", frame_count, ideal_gas + excess_gas, ideal_liquid + excess_liquid, gas_count, liquid_count);

            // reset averages for next round
            avg_count = 0;
            ideal_pot_gas_sum = 0.0;
            ideal_pot_liquid_sum = 0.0;
            widom_sum_gas = 0.0;
            widom_sum_liquid = 0.0;
        }

        if !trj_reader.update_with_next(&mut frame) {
            break;
        }
    }

}


fn get_particle_insertion_energy(rx: &[f64], ry: &[f64], rz: &[f64], num_particles: usize, x: f64, y: f64, z: f64, l_x: f64,l_y: f64, l_z: f64, cutoff_sqr: f64, e_shift: f64) -> f64 {
    let mut energy = 0.0;
    let half_l_x = l_x/2.0;
    let half_l_y = l_y/2.0;
    let half_l_z = l_z/2.0;
    for i in 0..num_particles {

        let dist_squared = get_particle_distance_squared(x, y, z, rx[i],ry[i],rz[i], l_x, l_y, l_z, half_l_x, half_l_y, half_l_z);
        if dist_squared < cutoff_sqr {
            energy += eval_pair_energy(dist_squared, e_shift).0;
        }
    }
    return energy;
}
