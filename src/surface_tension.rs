mod trajectory;
use trajectory::*;
mod energy;
use energy::*;
use std::env;

const LJ_EPS : f64 = 1.0;
const LJ_SIG : f64 = 1.0;

const AVG_OUTPUT_INTERVAL : usize = 10;

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
    println!("{:?}", frame);

    let box_half_x = frame.box_x / 2.0;
    let box_half_y = frame.box_y / 2.0;
    let box_half_z = frame.box_z / 2.0;


    let variable_without_name = frame.temperature/LJ_EPS * density;

    println!("Calculating surface tension");
    println!("~~~ THIS IS A RUNNING AVERAGE! ~~~");
    let mut trace_xy_sum = 0.0;
    let mut trace_z_sum = 0.0;
    let mut frame_count = 0;

    loop {
        frame_count += 1;

        let mut trace_xy = 0.0;
        let mut trace_z = 0.0;
        for i in 0..num_particles {
            for j in i+1..num_particles {
                // this needs some optimization for speed
                let dist_sqrt = get_particle_distance_squared(frame.rx[i], frame.ry[i], frame.rz[i], frame.rx[j], frame.ry[j], frame.rz[j], frame.box_x, frame.box_y, frame.box_z, box_half_x, box_half_y, box_half_z);
                let dist = dist_sqrt.sqrt();
                let dx = get_distance_with_pbc(frame.rx[i], frame.rx[j], frame.box_x, box_half_x);
                let dy = get_distance_with_pbc(frame.ry[i], frame.ry[j], frame.box_y, box_half_y);
                let dz = get_distance_with_pbc(frame.rz[i], frame.rz[j], frame.box_z, box_half_z);
                let virial = eval_virial(dist, LJ_EPS, LJ_SIG);
                trace_xy += (dx * dx + dy * dy) / dist * virial;
                trace_z += (dz * dz) / dist * virial;
            }
        }
        trace_xy_sum += trace_xy;
        trace_z_sum += trace_z;

        ///////////////////////////////////
        if frame_count % AVG_OUTPUT_INTERVAL == 0 {
            let p_z_avg = variable_without_name - 1.0/volume*(trace_z_sum/frame_count as f64);
            let p_xy_avg = variable_without_name - 1.0/2.0/volume*(trace_xy_sum/frame_count as f64);

            let surface_tension = eval_surface_tension(frame.box_z, p_z_avg, p_xy_avg);
            println!("Frame {}\t\tzz: {:.5}\txy: {:.5}\t\ttension: {:.5}", frame_count, p_z_avg, p_xy_avg, surface_tension);

        }

        // read next frame
        if !trj_reader.update_with_next(&mut frame) { break }
    }

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
