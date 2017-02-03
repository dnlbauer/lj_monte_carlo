
pub fn get_total_energy(rx: &[f64], ry: &[f64], rz: &[f64], num_particles: usize, l_x: f64, l_y: f64, l_z: f64, cutoff_squared: f64, e_corr: f64, e_shift: f64) -> (f64, f64) {
    let mut energy = 0.0;
    let mut virial = 0.0;
    let hl_x = l_x / 2.0;
    let hl_y = l_y / 2.0;
    let hl_z = l_z / 2.0;
    for i in 0..num_particles {
        for j in i+1..num_particles {
            let dist_squared = get_particle_distance_squared(rx[i], ry[i],rz[i],rx[j],ry[j],rz[j], l_x, l_y, l_z, hl_x, hl_y, hl_z);
            if dist_squared < cutoff_squared {
                let (e,v) = eval_pair_energy(dist_squared, e_shift);
                energy += e;
                virial += v;
            }
        }
    }
    energy += num_particles as f64 * e_corr;
    return (energy, virial);
}

pub fn get_particle_energy(rx: &[f64], ry: &[f64], rz: &[f64], p_index: usize, num_particles: usize, l_x: f64, l_y: f64, l_z: f64, cutoff_squared: f64, e_shift: f64) -> (f64, f64) {
    let mut energy = 0.0;
    let mut virial = 0.0;
    let hl_x = l_x / 2.0;
    let hl_y = l_y / 2.0;
    let hl_z = l_z / 2.0;
    for i in 0..num_particles {
        if i == p_index { continue; }

        let dist_squared = get_particle_distance_squared(rx[i], ry[i], rz[i], rx[p_index], ry[p_index], rz[p_index], l_x, l_y, l_z, hl_x, hl_y, hl_z);
        if dist_squared < cutoff_squared {
            let (e,v) = eval_pair_energy(dist_squared, e_shift);
            energy += e;
            virial += v;
        }
    }
    return (energy, virial);
}

pub fn get_particle_distance_squared(x1: f64,y1: f64,z1: f64,x2: f64,y2: f64,z2: f64, l_x: f64, l_y: f64, l_z: f64, hl_x: f64, hl_y: f64, hl_z: f64) -> f64 {
    let mut dx = (x1 - x2).abs();
    let mut dy = (y1 - y2).abs();
    let mut dz = (z1 - z2).abs();

    if dx > hl_x { dx -= l_x }
        else if dx < -hl_x { dx += l_x }
    if dy > hl_y { dy -= l_y}
        else if dy < -hl_y{ dy += l_y }
    if dz > hl_z { dz -= l_y }
        else if dz < -hl_z { dz += l_z}
    return dx*dx + dy*dy + dz*dz;

}

pub fn get_particle_distance(x1: f64,y1: f64,z1: f64,x2: f64,y2: f64,z2: f64, l_x: f64, l_y: f64, l_z: f64, hl_x: f64, hl_y: f64, hl_z: f64) -> f64{
    return get_particle_distance_squared(x1,y1,z1,x2,y2,z2, l_x, l_y, l_z, hl_x, hl_y, hl_z).sqrt();
}

#[test]
fn test_get_particle_distance_squared() {
    let (x1, y1, z1) = (0.0, 0.0, 0.0);
    let (x2, y2, z2) = (5.0, 0.0, 0.0);

    // no pbc
    assert!( (get_particle_distance_squared(x1,y1,z1,x2,y2,z2, 20.0, 20.0, 20.0, 10.0, 10.0, 10.0) - 25.0) < 0.00001);

    // with pbc
    assert!( (get_particle_distance_squared(x1,y1,z1,x2,y2,z2, 9.0,9.0,9.0, 4.5,4.5,4.5) - 16.0) < 0.00001, "{}", get_particle_distance_squared(x1,y1,z1,x2,y2,z2, 9.0,9.0,9.0, 4.5,4.5,4.5));

    // PBC edge case
    let (x3, y3, z3) = ( 0.948369102634727,1.4018028642626956,2.6884871697542323);
    let (x4, y4, z4) = ( 0.4924652252404308,0.27903672240586597,1.258265555104697);
    let a = get_particle_distance_squared(x3,y3,z3,x4,y4,z4, 1.418983411970384,1.418983411970384,2.83796682394076,0.709491705985192,0.709491705985192,1.418983411970384);
    let b = get_particle_distance_squared(x4,y4,z4,x3,y3,z3, 1.418983411970384,1.418983411970384,2.83796682394076,0.709491705985192,0.709491705985192,1.418983411970384);
    assert!( (a-b).abs() < 0.0000000001);

    let (x1, y1, z1) = (1.0, 1.0, 1.0);
    let (x2, y2, z2) = (99.0, 99.0, 99.0);
    let dist = get_particle_distance_squared(x1,y1,z1,x2,y2,z2, 100.0, 100.0, 100.0, 50.0, 50.0, 50.0);
    assert!(dist - 12.0 < 0.00001);
}

pub fn eval_pair_energy(dist_squared: f64, e_shift: f64) -> (f64, f64) {
    let r6 = ::LJ_SIG/(dist_squared * dist_squared * dist_squared);
    let r62 = r6*r6;
    let energy = 4.0 * ::LJ_EPS * (r62 - r6) - e_shift;
    let virial = 48.0 * ::LJ_EPS * ( r62 - 0.5 * r6 );
    return (energy, virial);
}

#[test]
fn test_eval_pair_energy() {
    let (e,v) = eval_pair_energy(1.0, 0.0);
    assert!( (e - 0.0).abs() < 0.00001, "{}",  e);

    let (e,v) = eval_pair_energy(2.0, 0.0);
    assert!( (e - -0.4375).abs() < 0.00001, "{}",  e);

    let (e,v) = eval_pair_energy(0.5, 0.0);
    assert!( (e - 224.0).abs() < 0.00001, "{}",  e);
}


pub fn eval_virial(distance: f64, LJ_EPS: f64, LJ_SIG: f64) -> f64 {
    let r7 = (LJ_SIG/distance).powi(7);
    let r13 = (LJ_SIG/distance).powi(13);
    return 24.0 * LJ_EPS / LJ_SIG * ( r7-2.0*r13 );
}

#[test]
fn test_eval_virial() {
    let dist = 1.5;
    let result = eval_virial(dist, 1.0, 1.0, 0.0);
    let expected = 1.1580288;
    assert!( (result - expected).abs() < 0.0001, "{}", result );
}
