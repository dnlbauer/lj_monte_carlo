mod trajectory;
use trajectory::*;
mod energy;
use energy::*;

const LJ_EPS : f64 = 1.0;
const LJ_SIG : f64 = 1.0;


fn get_virial(distance_sqr: f64) -> f64 {
    let r2 = LJ_SIG.powi(2)/distance_sqr;
    let r6 = r2 * r2 * r2;
    return 48.0 * LJ_EPS / LJ_SIG * ( r6*r6 - 0.5*r6 )
}

fn get_distance_with_pbc(x1: f64, x2: f64, length: f64, half_length: f64) -> f64 {
    let mut d = x1-x2;
    if d > half_length { d -= length }
    else if d < -half_length { d += length }
    return d;
}

fn main() {

    let mut trj_reader = TrjReader::new(&"2phases/1kk_100step.xyz".to_string());
    let mut frame = trj_reader.next_frame();

    let volume = frame.box_x * frame.box_y * frame.box_z;
    let density = frame.num_particles as f64 / volume;
    let num_particles = frame.num_particles;

    let box_half_x = frame.box_x / 2.0;
    let box_half_y = frame.box_y / 2.0;
    let box_half_z = frame.box_z / 2.0;

    let mut frame_count = 0;
    let mut p_xy_sum = 0.0;
//    let mut p_y_sum = 0.0;
    let mut p_z_sum = 0.0;

    let variable_without_name = frame.temperature/LJ_EPS * density;

    loop {
        frame_count += 1;

        let mut trace_xy = 0.0;
//        let mut trace_y = 0.0;
        let mut trace_z = 0.0;
        for i in 0..num_particles {
            for j in i+1..num_particles {
                let dist_sqrt = get_particle_distance_squared(frame.rx[i], frame.ry[i], frame.rz[i], frame.rx[j], frame.ry[j], frame.rz[j], frame.box_x, frame.box_y, frame.box_z, box_half_x, box_half_y, box_half_z);
                let dist = dist_sqrt.sqrt();
                let dx = get_distance_with_pbc(frame.rx[i], frame.rx[j], frame.box_x, box_half_x);
                let dy = get_distance_with_pbc(frame.ry[i], frame.ry[j], frame.box_y, box_half_y);
                let dz = get_distance_with_pbc(frame.rz[i], frame.rz[j], frame.box_z, box_half_z);
                let virial = get_virial(dist_sqrt);
                trace_xy += (dx * dx + dy * dy) / dist * virial;
//                trace_y += (dy * dy) / dist * virial;
                trace_z += (dz * dz) / dist * virial;
            }
        }
        let p_xy = variable_without_name - 1.0/(2.0*volume)*(trace_xy/num_particles as f64);
//        let p_yy = variable_without_name - 1.0/volume*(trace_y/num_particles);
        let p_zz = variable_without_name - 1.0/volume*(trace_z/num_particles as f64);

        p_xy_sum += p_xy;
//        p_y_sum += p_yy;
        p_z_sum += p_zz;

        ///////////////////////////////////

        if frame_count % 100 == 0 {
            print!(".");
        }

        if frame_count % 100 == 0 {
            let p_z_avg = p_z_sum / frame_count as f64;
            let p_xy_avg = p_xy_sum / frame_count as f64;
            let p_diff = p_z_avg - p_xy_avg;
            let surface_tension = (frame.box_z/2.0)*p_diff;
            println!("{}  zz: {}   xy: {}    diff: {}    tension:   {}", frame_count, p_z_avg, p_xy_avg, p_diff, surface_tension);

            frame_count = 0;
            p_z_sum = 0.0;
            p_xy_sum = 0.0;
//            trace_xy_sum = 0.0;
//            trace_z_sum = 0.0;
//            p_tangial_sum = 0.0;
//            p_normal_sum = 0.0

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