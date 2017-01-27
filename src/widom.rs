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
const RUN_AVG_SIZE : usize = 10;

const SHIFT : bool = true;

fn main() {
    // open file
    let args: Vec<String> = env::args().collect();
    let mut trj_reader = TrjReader::new(&args[1]);

    let mut frame = trj_reader.next_frame();
    println!("{:?}", frame);
    // get some non changing values
    let volume = frame.box_x * frame.box_y * frame.box_z;
    let beta = 1.0/frame.temperature;
    let cutoff_sqr = frame.lj_cutoff * frame.lj_cutoff;

    // account for gas phase slab
    let gas_slab = args[2].parse::<f64>().unwrap();
    let liquid_height = frame.box_z / (gas_slab + 1.0);
    let liquid_volume = volume / (gas_slab + 1.0);
    let gas_volume = volume - liquid_volume;

    let wave = MKSA_PLANCKS_CONSTANT_H / (2.0 * std::f64::consts::PI * MASS * frame.temperature / frame.lj_eps).sqrt();
    let e_shift = if SHIFT { 4.0 * LJ_EPS * ( (LJ_SIG/frame.lj_cutoff).powi(12) - (LJ_SIG/frame.lj_cutoff).powi(6) ) } else { 0.0 };

    let mut rng = rand::thread_rng();


    // average counters
    let mut frame_count = 0;
    let mut widom_sum_gas = 0.0;
    let mut ideal_sum_gas = 0.0;
    let mut widom_sum_liquid = 0.0;
    let mut ideal_sum_liquid = 0.0;

    // loop over frames
    let mut avg_count = 0;
    loop  {
        frame_count += 1;

        // calculate number of particles in each phase
        let mut liquid_count = 0.0;
        let mut gas_count = 0.0;
        for i in 0..frame.num_particles {
            if frame.rz[i] < liquid_height { liquid_count += 1.0; }
            else { gas_count += 1.0; }
        }

        // test particle insertion multiple times
        for i in 0..INSERTIONS_PER_STEP {
            avg_count += 1;

            // liquid test partciles
            let lx = frame.box_x * rng.gen::<f64>();
            let ly = frame.box_y * rng.gen::<f64>();
            let lz = liquid_height * rng.gen::<f64>();

            let widom_e_liquid = get_particle_insertion_energy(&frame.rx, &frame.ry, &frame.rz, frame.num_particles, lx, ly, lz, frame.box_x, frame.box_y, frame.box_z, cutoff_sqr, e_shift);
            widom_sum_liquid += (-beta*widom_e_liquid).exp();

            // test particle gas energy
            let gx = frame.box_x * rng.gen::<f64>();
            let gy = frame.box_y * rng.gen::<f64>();
            let gz = ((frame.box_z - liquid_height) * rng.gen::<f64>()) + liquid_height;
            let widom_e_gas = get_particle_insertion_energy(&frame.rx, &frame.ry, &frame.rz, frame.num_particles, gx, gy, gz, frame.box_x, frame.box_y, frame.box_z, cutoff_sqr, e_shift);
            widom_sum_gas += (-beta*widom_e_gas).exp();

            // calculate ideal gas potential for both phases
            ideal_sum_gas += frame.temperature / frame.lj_eps * (gas_volume/(wave.powi(3)* gas_count)).ln();
            ideal_sum_liquid += frame.temperature / frame.lj_eps * (liquid_volume/(wave.powi(3)* liquid_count)).ln();

        }

        // print running averages
        if avg_count / INSERTIONS_PER_STEP > RUN_AVG_SIZE {
            let ideal_gas = ideal_sum_gas / avg_count as f64;
            let ideal_liquid = ideal_sum_liquid / avg_count as f64;
            let excess_gas = -(widom_sum_gas/avg_count as f64).ln()/beta;
            let excess_liquid = -(widom_sum_liquid/avg_count as f64).ln()/beta;
            println!("Frame {}\tgas {}\tliquid {}\t particles gas/liquid:{}/{}", frame_count, ideal_gas + excess_gas, ideal_liquid + excess_liquid, gas_count, liquid_count);

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
