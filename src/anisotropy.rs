mod trajectory;
use trajectory::*;
mod energy;
use energy::*;
use std::env;

const LJ_EPS : f64 = 1.0;
const LJ_SIG : f64 = 1.0;

const AVG_OUTPUT_INTERVAL : usize = 10;
const SLAB_NUM : usize = 200;

fn main() {
    let mut filename = "montecarlo.xyz".to_string();
    let mut skip: usize = 0;

    // parse cmd line args
    let args: Vec<String> = env::args().collect();
    for i in 0..args.len() {
        if args[i] == "-f" {
            filename = args[i + 1].clone();
        } else if args[i] == "-s" {
            skip = args[i + 1].parse::<usize>().unwrap();
        }
    }

    // open file and skip to requested position
    let mut trj_reader = TrjReader::new(&filename);
    if skip > 0 { trj_reader.skip(skip) };

    // read first trajectory and system params
    let mut frame = trj_reader.next_frame();
    let volume = frame.box_x * frame.box_y * frame.box_z;
    let density = frame.num_particles as f64 / volume;
    let num_particles = frame.num_particles;
    println!("#{:?}", frame);

    let box_half_x = frame.box_x / 2.0;
    let box_half_y = frame.box_y / 2.0;
    let box_half_z = frame.box_z / 2.0;

    let mut frame_count = 0;
    let mut p_xy_sum = 0.0;
    let mut p_z_sum = 0.0;

    let slab_height = frame.box_z / SLAB_NUM as f64;
    let slab_volume = frame.box_x * frame.box_y * slab_height;
    let mut virial_histogram_xy = [0.0; SLAB_NUM];
    let mut virial_histogram_xy_counter = [0; SLAB_NUM];
    let mut virial_histogram_z = [0.0; SLAB_NUM];
    let mut virial_histogram_z_counter = [0; SLAB_NUM];
    let mut slab_particles_sum = [0; SLAB_NUM];

    loop {
        frame_count += 1;

        for i in 0..num_particles {
            let slab_i = get_slab_number_for_position(frame.rz[i], slab_height);
            slab_particles_sum[slab_i-1] += 1;
            for j in i+1..num_particles {
                // distance
                let dist_sqrt = get_particle_distance_squared(frame.rx[i], frame.ry[i], frame.rz[i], frame.rx[j], frame.ry[j], frame.rz[j], frame.box_x, frame.box_y, frame.box_z, box_half_x, box_half_y, box_half_z);
                let dist = dist_sqrt.sqrt();
                let dx = get_distance_with_pbc(frame.rx[i], frame.rx[j], frame.box_x, box_half_x);
                let dy = get_distance_with_pbc(frame.ry[i], frame.ry[j], frame.box_y, box_half_y);
                let dz = get_distance_with_pbc(frame.rz[i], frame.rz[j], frame.box_z, box_half_z);

                // calculate virial tensor
                let virial = eval_virial(dist, LJ_EPS, LJ_SIG);
                let virial_xy = (dx * dx + dy * dy) / dist.powi(2) * virial;
                let virial_z = (dz * dz) / dist * virial;

                // calculate slab distribution
                let slab_j = get_slab_number_for_position(frame.rz[j], slab_height);
                let first_slab_index = get_first_slab_for_trace(frame.rz[i], frame.rz[j], slab_height) - 1;
                let last_slab_index = get_last_slab_for_trace(frame.rz[i], frame.rz[j], slab_height) - 1 ;
                let num_slabs = last_slab_index - first_slab_index + 1;
                // println!("{} {} {}", first_slab_index, last_slab_index, num_slabs);
                let virial_xy_partial = virial_xy / num_slabs as f64;
                let virial_z_partial = virial_z / num_slabs as f64;

                for slab in first_slab_index..last_slab_index+1 {
                    virial_histogram_xy[slab] += virial_xy_partial;
                    virial_histogram_xy_counter[slab] += 1;
                    virial_histogram_z[slab] += virial_z_partial;
                    virial_histogram_z_counter[slab] += 1;
                }
            }
        }
        println!("# Frame {}", frame_count);
        // produce output
        if frame_count > 10 {
            println!("# Slab\tdensity\txy\tz\tanisotropy");
            for i in 0..SLAB_NUM {
                let slab_density = slab_particles_sum[i] as f64 / frame_count as f64 / slab_volume;
                let variable_without_name = frame.temperature/LJ_EPS * slab_density;
                let p_xy = variable_without_name - 1.0/(2.0*volume)*( virial_histogram_z[i]);
                let p_zz = variable_without_name - 1.0/volume*( virial_histogram_xy[i]);
                let anisotropy = p_zz - p_xy;
                println!("{} {} {} {} {}", i+1, slab_density, p_xy, p_zz, anisotropy);
            }
            std::process::exit(0);
        }
        //

        //
        // p_xy_sum += p_xy;
        // p_z_sum += p_zz;
        //
        // ///////////////////////////////////
        // if frame_count % AVG_OUTPUT_INTERVAL == 0 {
        //     let p_z_avg = p_z_sum / frame_count as f64;
        //     let p_xy_avg = p_xy_sum / frame_count as f64;
        //     let p_diff = p_z_avg - p_xy_avg;
        //     let surface_tension = eval_surface_tension(frame.box_z, p_z_avg, p_xy_avg);
        //     println!("Frame {}\t\tzz: {:.5}\txy: {:.5}\tdifference: {:.5}\t\ttension: {:.5}", frame_count, p_z_avg, p_xy_avg, p_diff, surface_tension);
        //
        //     // frame_count = 0;
        //     // p_z_sum = 0.0;
        //     // p_xy_sum = 0.0;
        // }

        // read next frame
        if !trj_reader.update_with_next(&mut frame) { break }
    }

}

fn get_first_slab_for_trace(rz1: f64, rz2:f64, slab_height: f64) -> usize {
    let first = rz1.min(rz2);
    let mut slab_num = 1;
    let mut slab_max_height = slab_height;
    while(slab_max_height <= first) {
        slab_num += 1;
        slab_max_height += slab_height;
    }
    return slab_num;
}


#[test]
fn test_get_first_slab_for_trace() {
    let rz1 = 2.5;
    let rz2 = 3.7;
    let slab_height = 1.1;
    let expected = 3;
    let result = get_first_slab_for_trace(rz1, rz2, slab_height);
    assert_eq!(expected, result, "{}", result);

    let rz1 = 3.7;
    let rz2 = 2.5;
    let slab_height = 1.1;
    let expected = 3;
    let result = get_first_slab_for_trace(rz1, rz2, slab_height);
    assert_eq!(expected, result, "{}", result);

    let rz1 = 0.0;
    let rz2 = 2.5;
    let slab_height = 1.1;
    let expected = 1;
    let result = get_first_slab_for_trace(rz1, rz2, slab_height);
    assert_eq!(expected, result, "{}", result);

    let rz1 = 1.1;
    let rz2 = 0.0;
    let slab_height = 1.1;
    let expected = 1;
    let result = get_first_slab_for_trace(rz1, rz2, slab_height);
    assert_eq!(expected, result, "{}", result);

    let rz1 = 1.0;
    let rz2 = 10.0;
    let slab_height = 1.0;
    let expected = 2;
    let result = get_first_slab_for_trace(rz1, rz2, slab_height);
    assert_eq!(expected, result, "{}", result);
}

fn get_last_slab_for_trace(rz1: f64, rz2:f64, slab_height: f64) -> usize {
    let last = rz1.max(rz2);
    let mut slab_num = 1;
    let mut slab_max_height = slab_height;
    while(slab_max_height < last) {
        slab_num += 1;
        slab_max_height += slab_height;
    }
    return slab_num;

}

#[test]
fn test_get_last_slab_for_trace() {
    let rz1 = 2.5;
    let rz2 = 3.7;
    let slab_height = 1.1;
    let expected = 4;
    let result = get_last_slab_for_trace(rz1, rz2, slab_height);
    assert_eq!(expected, result, "{}", result);

    let rz1 = 3.7;
    let rz2 = 2.5;
    let slab_height = 1.1;
    let expected = 4;
    let result = get_last_slab_for_trace(rz1, rz2, slab_height);
    assert_eq!(expected, result, "{}", result);

    let rz1 = 0.0;
    let rz2 = 2.5;
    let slab_height = 1.1;
    let expected = 3;
    let result = get_last_slab_for_trace(rz1, rz2, slab_height);
    assert_eq!(expected, result, "{}", result);

    let rz1 = 1.1;
    let rz2 = 0.0;
    let slab_height = 1.1;
    let expected = 1;
    let result = get_last_slab_for_trace(rz1, rz2, slab_height);
    assert_eq!(expected, result, "{}", result);

    let rz1 = 1.0;
    let rz2 = 10.0;
    let slab_height = 1.0;
    let expected = 10;
    let result = get_last_slab_for_trace(rz1, rz2, slab_height);
    assert_eq!(expected, result, "{}", result);
}

pub fn get_slab_number_for_position(z: f64, slab_height: f64) -> usize {
    let mut slab = 1;
    let mut z_counter = slab_height;
    while z_counter <= z {
        z_counter += slab_height;
        slab += 1;
    }
    return slab;
}

#[test]
fn test_get_slab_no() {
    let z = 3.5;
    let slab_height = 1.0;
    assert_eq!(4, get_slab_number_for_position(z, slab_height));

    let z = 0.0;
    let slab_height = 1.0;
    assert_eq!(1, get_slab_number_for_position(z, slab_height));

    let z = 0.1;
    let slab_height = 0.2;
    assert_eq!(1, get_slab_number_for_position(z, slab_height));

    let z = 0.0001;
    let slab_height = 0.2;
    assert_eq!(1, get_slab_number_for_position(z, slab_height));


    let z = 0.19999999_f64;
    let slab_height = 0.2;
    assert_eq!(1, get_slab_number_for_position(z, slab_height));

    let z = 1.0;
    let slab_height = 1.0;
    assert_eq!(2, get_slab_number_for_position(z, slab_height));
}


/// calc surface tension from box z size and pressure tensor
fn eval_surface_tension(box_z: f64, p_zz: f64, p_xy: f64) -> f64 {
    return box_z / 2.0 * (p_zz - p_xy);
}

#[test]
fn test_eval_surface_tension() {
    let expected = 2.0;
    let result = eval_surface_tension(2.0,5.0,3.0);
    assert!( (result-expected).abs() < 0.0001, "{}", result );
}
