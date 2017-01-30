#![allow(dead_code)]
#![allow(unused_must_use)] // hate.
#![allow(unused_variables)]

use std::error::Error;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use std::fmt;
use std::io::BufReader;


pub struct XYZTrajectory {
    file: File,

}

impl XYZTrajectory {
    pub fn new(filename: &String) -> XYZTrajectory {
        let path = Path::new(filename);
        let display = path.display();
        // Open a file in write-only mode, returns `io::Result<File>`
        let traj_file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}",
                               display,
                               why.description()),
            Ok(file) => file,
        };

        XYZTrajectory { file: traj_file }

    }


    pub fn write(&mut self, rx: &[f64], ry: &[f64], rz: &[f64], num_particles: usize, box_x : f64, box_y : f64, box_z : f64, temp: f64 ,lj_eps : f64, lj_sig : f64, lj_cutoff : f64, flush: bool) {
        self.file.write(format!("{} ## Box: {} {} {} Temp: {} LJ: {}/{}/{}\n", num_particles, box_x,box_y,box_z,temp, lj_eps, lj_sig, lj_cutoff).as_bytes());
        for i  in 0..num_particles {
            let formatted = format!("atom{} {} {} {}\n",
                                    i+1, rx[i], ry[i], rz[i]
            );
            self.file.write(formatted.as_bytes());
        }

        if flush { self.file.flush(); }

    }
}

pub struct Frame {
    pub rx : Vec<f64>,
    pub ry : Vec<f64>,
    pub rz : Vec<f64>,
    pub num_particles: usize,
    pub box_x: f64,
    pub box_y: f64,
    pub box_z: f64,
    pub temperature: f64,
    pub lj_eps: f64,
    pub lj_sig: f64,
    pub lj_cutoff: f64,
}

impl fmt::Debug for Frame {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("Frame")
            .field("num_particles", &self.num_particles)
            .field("box_x", &self.box_x)
            .field("box_y", &self.box_y)
            .field("box_z", &self.box_z)
            .field("temperature", &self.temperature)
            .field("lj_eps", &self.lj_eps)
            .field("lj_sig", &self.lj_sig)
            .field("lj_cutoff", &self.lj_cutoff)
            .finish()
    }
}

pub struct TrjReader {
    pub reader: BufReader<File>,
}
impl TrjReader {
    pub fn new(filename: &String) -> TrjReader {
        let file = File::open(filename).expect("Failed to open file.");
        let reader : BufReader<File> = BufReader::new(file);
        return TrjReader { reader:reader };
    }

    // get next frame
    pub fn next_frame(&mut self) -> Frame {
        let buffer_string = &mut String::new();
        match self.reader.read_line(buffer_string) {
            Ok(size) => {
                if size == 0 {
                    panic!("no new frame!")
                }
            },
            Err(e) => panic!("no new frame!")
        }
        let first_line_vec : Vec<&str> = buffer_string.split_whitespace().collect();
        let num_particles = first_line_vec[0].parse::<usize>().unwrap();
        let box_x = first_line_vec[3].parse::<f64>().unwrap();
        let box_y = first_line_vec[4].parse::<f64>().unwrap();
        let box_z = first_line_vec[5].parse::<f64>().unwrap();
        let temp = first_line_vec[7].parse::<f64>().unwrap();
        let lj : Vec<&str> = first_line_vec[9].split("/").collect();
        let lj_eps = lj[0].parse::<f64>().unwrap();
        let lj_sig = lj[1].parse::<f64>().unwrap();
        let lj_cutoff = lj[2].parse::<f64>().unwrap();
        let mut rx = Vec::new();
        let mut ry = Vec::new();
        let mut rz = Vec::new();
        for i in 0..num_particles{
            let atom_line = &mut String::new();
            match self.reader.read_line(atom_line) {
                Ok(size) => {
                    let atom_vec : Vec<&str> = atom_line.split_whitespace().collect();
                    rx.push(atom_vec[1].parse::<f64>().unwrap());
                    ry.push(atom_vec[2].parse::<f64>().unwrap());
                    rz.push(atom_vec[3].parse::<f64>().unwrap());
                },
                Err(e) => return panic!("no new frame!")
            };
        }
        let frame = Frame {
            rx : rx,
            ry : ry,
            rz : rz,
            num_particles : num_particles,
            box_x : box_x,
            box_y : box_y,
            box_z : box_z,
            temperature : temp,
            lj_eps : lj_eps,
            lj_sig : lj_sig,
            lj_cutoff : lj_cutoff,
        };
        return frame;
    }

    // read next frame data into the frame
    pub fn update_with_next(&mut self, frame: &mut Frame) -> bool {
        let buffer_string = &mut String::new();
        match self.reader.read_line(buffer_string) {
            Ok(size) => {
                if size == 0 {
                    return false;
                }
            },
            Err(e) => { return false; }
        }
        let first_line_vec : Vec<&str> = buffer_string.split_whitespace().collect();
        frame.num_particles = first_line_vec[0].parse::<usize>().unwrap();
        frame.box_x = first_line_vec[3].parse::<f64>().unwrap();
        frame.box_y = first_line_vec[4].parse::<f64>().unwrap();
        frame.box_z = first_line_vec[5].parse::<f64>().unwrap();
        frame.temperature = first_line_vec[7].parse::<f64>().unwrap();
        let lj : Vec<&str> = first_line_vec[9].split("/").collect();
        frame.lj_eps = lj[0].parse::<f64>().unwrap();
        frame.lj_sig = lj[1].parse::<f64>().unwrap();
        frame.lj_cutoff = lj[2].parse::<f64>().unwrap();
        for i in 0..frame.num_particles {
            let atom_line = &mut String::new();
            match self.reader.read_line(atom_line) {
                Ok(size) => {
                    if size == 0 {
                        return false;
                    } else {
                        let atom_vec : Vec<&str> = atom_line.split_whitespace().collect();
                        frame.rx[i] = atom_vec[1].parse::<f64>().unwrap();
                        frame.ry[i]= atom_vec[2].parse::<f64>().unwrap();
                        frame.rz[i] = atom_vec[3].parse::<f64>().unwrap();
                    }
                },
                Err(e) => {return false;}
            }
        }

        return true;
    }

    // skip x frames
    pub fn skip(&mut self, skip: usize) {
        if skip < 1 { return };

        // find number of particles
        let buffer_string = &mut String::new();
        match self.reader.read_line(buffer_string) {
            Ok(size) => {
                if size == 0 {
                    panic!("no new frame!")
                }
            },
            Err(e) => panic!("no new frame!")
        }
        let first_line_vec : Vec<&str> = buffer_string.split_whitespace().collect();
        let num_particles = first_line_vec[0].parse::<usize>().unwrap();

        let lines_to_skip = (num_particles + 1) * skip - 1;
        let mut skipped = 0;
        loop {
            self.reader.read_line(&mut String::new());
            skipped += 1;
            if skipped == lines_to_skip { break; }
        }
    }

}
