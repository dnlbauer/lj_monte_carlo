mod trajectory;
use trajectory::*;
mod energy;
use energy::*;

const LJ_EPS : f64 = 1.0;
const LJ_SIG : f64 = 1.0;


fn get_derived_pair_potential(distance: f64) -> f64 {
    let r = LJ_SIG/distance;
    let r3 = r*r*r;
    let r4 = r*r*r*r;
    return 24.0 * LJ_EPS / LJ_SIG * ( r4*r3 - 2.0*(r4*r4*r4*r) )
}

fn get_distance_with_pbc(x1: f64, x2: f64, length: f64, half_length: f64) -> f64 {
    let mut d = x1-x2;
    if d > half_length { d -= length }
    else if d < -half_length { d += length }
    return d;
}

fn main() {

    let mut trj_reader = TrjReader::new(&"2phases/10k_1step.xyz".to_string());
    let mut frame = trj_reader.next_frame();

    let volume = frame.box_x * frame.box_y * frame.box_z;
    let density = frame.num_particles as f64 / volume;

    let box_half_x = frame.box_x / 2.0;
    let box_half_y = frame.box_y / 2.0;
    let box_half_z = frame.box_z / 2.0;

    let mut frame_count = 0;
    let mut trace_xy_sum = 0.0;
    let mut trace_z_sum = 0.0;
    let mut p_normal_sum = 0.0;
    let mut p_tangial_sum = 0.0;
    let mut p_diff_sum = 0.0;

    let variable_without_name = frame.temperature/LJ_EPS * density;
    loop {
        frame_count += 1;

        let mut trace_xy = 0.0;
        let mut trace_z = 0.0;
        for i in 0..frame.num_particles {
            for j in i+1..frame.num_particles {
                let dist = get_particle_distance(frame.rx[i], frame.ry[i], frame.rz[i], frame.rx[j], frame.ry[j], frame.rz[j], frame.box_x, frame.box_y, frame.box_z, box_half_x, box_half_y, box_half_z);
                let dx = get_distance_with_pbc(frame.rx[i], frame.rx[j], frame.box_x, box_half_x);
                let dy = get_distance_with_pbc(frame.ry[i], frame.ry[j], frame.box_y, box_half_y);
                let dz = get_distance_with_pbc(frame.rz[i], frame.rz[j], frame.box_z, box_half_z);
                trace_xy += (dx * dx  + dy * dy) / dist * get_derived_pair_potential(dist);
                trace_z += (dz * dz) / dist * get_derived_pair_potential(dist);
            }
        }
        let p_tangial = variable_without_name - 1.0/(2.0*volume)*trace_xy;
        let p_normal = variable_without_name - 1.0/volume*trace_z;
        let p_diff = p_normal - p_tangial_sum;
        p_diff_sum += p_diff;
        p_tangial_sum += p_tangial;
        p_normal_sum += p_normal;
        trace_xy_sum += trace_xy;
        trace_z_sum += trace_z;

        ///////////////////////////////////

        if frame_count % 100 == 0 {
            print!(".");
        }

        if frame_count % 100 == 0 {
//            let trace_xy_avg = trace_xy_sum / frame_count as f64;
//            let trace_z_avg = trace_z_sum / frame_count as f64;
            let p_tangial = p_tangial_sum / frame_count as f64;
            let p_normal = p_normal_sum / frame_count as f64;
            let surface_tension = frame.box_z/2.0*(p_diff_sum/frame_count as f64);
            println!("{}  tangial: {}   normal: {}    diff: {}    tension:   {}", frame_count, p_tangial, p_normal, p_diff_sum/frame_count as f64, surface_tension);

            frame_count = 0;
            trace_xy_sum = 0.0;
            trace_z_sum = 0.0;
            p_tangial_sum = 0.0;
            p_normal_sum = 0.0

        }

        if !trj_reader.update_with_next(&mut frame) { break }
    }

}

//        for slab in 0..NUM_SLABS {
//
//            // calculate slab density
//            let slab_end_z = slab_height * slab as f64;
//            let slab_start_z = slab_end_z - slab_height;
//            let mut slab_particle_count = 0;
//            for i in 0..frame.num_particles {
//                if slab_start_z > frame.rz[i] && frame.rz[i] < slab_end_z {
//                    slab_particle_count += 1;
//                }
//            }
//            let slab_density = slab_particle_count as f64 / slab_volume;
//
//            for i in 0..frame.num_particles {
//                for j in i+1..frame.num_particles {
//
//                }
//            }
//
//            println!("Slab {}\tDensity: {}", slab, slab_density);
//        }