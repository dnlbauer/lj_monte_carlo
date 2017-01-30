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

const INSERTIONS_PER_STEP : usize = 1000;
const RUN_AVG_SIZE : usize = 100;

const SHIFT : bool = true;



fn main() {
    // parse args
    let args: Vec<String> = env::args().collect();
    let mut filename = "montecarlo.xyz".to_string();
    let mut skip: usize = 0;

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
        }
    }

    // open file and skip to requiested position
    let mut trj_reader = TrjReader::new(&filename);
    if skip > 0 { trj_reader.skip(skip) };

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

    // lambda = h/sqrt(2*pi*m*kb*T) = (2*pi*m*kb*T/h^2)^3/2
    // this is lambda^3
    let wave3 =  ((2.0 * std::f64::consts::PI * MASS * frame.temperature / frame.lj_eps)/MKSA_PLANCKS_CONSTANT_H.powi(2)).powf(3.0/2.0);
    let e_shift = if SHIFT { 4.0 * LJ_EPS * ( (LJ_SIG/frame.lj_cutoff).powi(12) - (LJ_SIG/frame.lj_cutoff).powi(6) ) } else { 0.0 };


    // average counters
    let mut frame_count = 0;
    let mut widom_sum_gas = 0.0;
    let mut ideal_sum_gas = 0.0;
    let mut widom_sum_liquid = 0.0;
    let mut ideal_sum_liquid = 0.0;
    let mut avg_count = 0;
    loop  {
        frame_count += 1;

        // calculate number of particles in each phase
        let mut liquid_count = 0.0;
        let mut gas_count = 0.0;
        for i in 0..frame.num_particles {
            if frame.rz[i] > liquid_start && frame.rz[i] < liquid_end { liquid_count += 1.0; }
            else { gas_count += 1.0; }
        }


        // test particle insertion multiple times
        for i in 0..INSERTIONS_PER_STEP {
            avg_count += 1;

            // liquid test partcile
            let lx = frame.box_x * rng.gen::<f64>();
            let ly = frame.box_y * rng.gen::<f64>();
            let lz = liquid_start + (liquid_height * rng.gen::<f64>());
            let widom_e_liquid = get_particle_insertion_energy(&frame.rx, &frame.ry, &frame.rz, frame.num_particles, lx, ly, lz, frame.box_x, frame.box_y, frame.box_z, cutoff_sqr, e_shift);
            widom_sum_liquid += (-beta*widom_e_liquid).exp();
//            widom_sum_liquid += widom_e_liquid;

            // gas test particle
            let upper_or_lower = rng.gen::<bool>();
            let gx = frame.box_x * rng.gen::<f64>();
            let gy = frame.box_y * rng.gen::<f64>();
            let gz;
            if upper_or_lower {
                gz = rng.gen::<f64>() * liquid_start;
            } else {
                gz = rng.gen::<f64>() * (frame.box_z - liquid_end) + liquid_end;
            }
            let widom_e_gas = get_particle_insertion_energy(&frame.rx, &frame.ry, &frame.rz, frame.num_particles, gx, gy, gz, frame.box_x, frame.box_y, frame.box_z, cutoff_sqr, e_shift);
            widom_sum_gas += (-beta*widom_e_gas).exp();
//            widom_sum_gas += widom_e_gas;

            // calculate ideal gas potential
            ideal_sum_gas += -frame.temperature / frame.lj_eps * ((gas_volume/gas_count) * wave3).ln();
            ideal_sum_liquid += frame.temperature / frame.lj_eps * ((liquid_volume/liquid_count) * wave3).ln();

        }

        // print running averages
        if avg_count / INSERTIONS_PER_STEP > RUN_AVG_SIZE {
            let ideal_gas = ideal_sum_gas / avg_count as f64;
            let ideal_liquid = ideal_sum_liquid / avg_count as f64;
//            let excess_gas = - (widom_sum_gas/avg_count as f64).ln()/beta;
//            let excess_liquid = - (widom_sum_liquid/avg_count as f64).ln()/beta;
            let excess_gas = -(widom_sum_gas/avg_count as f64).ln()/beta;
            let excess_liquid = -(widom_sum_liquid/avg_count as f64).ln()/beta;
            let gas_total = ideal_gas + excess_gas;
            let liquid_total = ideal_liquid + excess_liquid;
            println!("Frame {}\tg_ex: {:5}\tl_ex: {:5}\tg_tot: {:5}\tl_tot: {:5}\t\tnparticles: {}/{}", frame_count, excess_gas, excess_liquid, gas_total, liquid_total, gas_count, liquid_count);
//            println!("Frame {}\tgas {}\tliquid {}\t particles gas/liquid:{}/{}", frame_count, ideal_gas + excess_gas, ideal_liquid + excess_liquid, gas_count, liquid_count);

            // reset averages for next round
            avg_count = 0;
            ideal_sum_gas = 0.0;
            ideal_sum_liquid = 0.0;
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
