mod trajectory;
use trajectory::*;
use std::env;

fn main() {
    // open file
    let args: Vec<String> = env::args().collect();

    let mut filename : String = "montecarlo.xyz".to_string();
    let mut skip_frames : usize = 0;
    let mut slabs : usize = 256;
    for i in 0..args.len() {
        if args[i] == "-f" {
            filename = args[i+1].clone();
        } else if args[i] == "-s" {
            skip_frames = args[i+1].parse::<usize>().unwrap();
        } else if args[i] == "--slabs" {
            slabs = args[i+1].parse::<usize>().unwrap();
        }
    }

    let mut trj_reader = TrjReader::new(&filename);

    // skip some frames
    println!("# Skipping {} frames.", skip_frames);
    trj_reader.skip(skip_frames);
    let mut frame = trj_reader.next_frame();

    let slab_height = frame.box_z / slabs as f64;
    let slab_volume = slab_height * frame.box_x * frame.box_y;
    println!("# Density calculation with {} slabs (height={})", slabs, slab_height);

    let mut slab_particles : Vec<f64> = vec![0.0; slabs];
    let mut frame_count : usize = 0;

    // loop over frames
    loop  {
        frame_count += 1;
        for i in 0..frame.num_particles {
            let slab_no : usize = get_slab_number_for_position(frame.rz[i], slab_height) - 1;
            slab_particles[slab_no] += 1.0;
        }

        if !trj_reader.update_with_next(&mut frame) {
            break;
        }
    }

    println!("# Averaged over {} frames", frame_count);
    println!("# Position    Density");
    for i in 0..slab_particles.len() {
        let p_max = slab_height * (i as f64 + 1.0);
        let position =(p_max + p_max-slab_height) / 2.0;
        let particles = slab_particles[i] / frame_count as f64;
        let density = particles / slab_volume;
        println!("{}\t{}\t{}", position, density, particles);
    }

}

pub fn get_slab_number_for_position(z: f64, slab_height: f64) -> usize {
    let mut slab = 0;
    let mut z_counter = 0.0_f64;
    while z_counter  <= z {
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