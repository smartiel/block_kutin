use super::matrix::MatrixF2;
use super::topology::Topology;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};

fn reconstruct_genealogy(database: &HashMap<u64, (i32, u64)>, repr: &u64) -> Vec<u64> {
    let mut genealogy = Vec::new();
    let mut current = *repr;
    loop {
        // println!("{current} has parent :");
        let (_, parent) = database.get(&current).unwrap();
        // println!("{parent}");
        if *parent == 0 {
            break;
        }
        genealogy.push(*parent);
        current = *parent;
    }
    return genealogy;
}

/// Generates the circuit implementing a given matrix
/// /!\ This function modifies the matrix /!\
pub fn get_circuit(
    target_db: &HashMap<u64, Vec<u64>>,
    all_operators: &Vec<Vec<(usize, usize)>>,
    matrix: &mut MatrixF2,
) -> Vec<(usize, usize)> {
    let repr = matrix.cannonical_representation();
    let mut circuit = Vec::new();
    let history = target_db.get(&repr).expect("The operator is not valid!");
    for ancestor in history {
        for operator in all_operators.iter() {
            matrix.apply_operator(operator);
            let new_repr = matrix.cannonical_representation();
            if *ancestor == new_repr {
                circuit.extend_from_slice(operator);
                break;
            }
            matrix.apply_operator(operator);
        }
    }
    return circuit;
}

pub fn get_circuit_from_repr(
    target_db: &HashMap<u64, Vec<u64>>,
    all_operators: &Vec<Vec<(usize, usize)>>,
    repr: &u64,
    n: usize,
    m: usize,
) -> Vec<(usize, usize)> {
    let mut matrix = MatrixF2::from_representation(&n, &m, repr);
    return get_circuit(target_db, all_operators, &mut matrix);
}

pub fn enumerate_step_2(topology: &Topology) -> HashMap<u64, Vec<u64>> {
    let mut global_database: HashMap<u64, (i32, u64)> = HashMap::new();
    let mut target_database: HashMap<u64, Vec<u64>> = HashMap::new();
    let all_operators: Vec<Vec<(usize, usize)>> = topology.get_all_operators();
    let mut initial = MatrixF2::identity(&topology.nvertices, &topology.nvertices);
    let mut depth = 0;
    global_database.insert(initial.cannonical_representation(), (0, 0));
    let mut queue: Vec<(MatrixF2, u64)> = vec![(initial, 0)];
    // This is the total number of target operators that exist
    let target_count =
        (2 as usize).pow((topology.nvertices / 2) as u32 * (topology.nvertices / 2) as u32);
    // This will count the number of observed target operators (we stop when it reaches `target_count`)
    let mut count: usize = 0;
    loop {
        depth += 1;
        if queue.is_empty() {
            // This happens if we observed all possible operators (if should never trigger)
            break;
        }
        if count == target_count {
            // This happens if we observed every possible target operator
            break;
        }
        let mut exit = false;
        let mut new_queue: Vec<(MatrixF2, u64)> = Vec::new();
        for (mut matrix, _) in queue.drain(..) {
            let new_parent = matrix.cannonical_representation();
            for operator in all_operators.iter() {
                matrix.apply_operator(operator);
                let repr = matrix.cannonical_representation();
                if !global_database.contains_key(&repr) {
                    global_database.insert(repr, (depth, new_parent));
                    new_queue.push((matrix.clone(), new_parent));
                    if matrix.check_step_2() {
                        let genealogy = reconstruct_genealogy(&global_database, &repr);
                        target_database.insert(repr, genealogy);
                        count += 1;
                        if target_count == count {
                            exit = true;
                            break;
                        }
                    }
                }
                matrix.apply_operator(operator);
            }
            if exit {
                break;
            }
        }
        queue = new_queue;
    }

    return target_database;
}

pub fn enumerate_step_1(topology: &Topology) -> HashMap<u64, Vec<u64>> {
    let mut global_database: HashMap<u64, (i32, u64)> = HashMap::new();
    let mut target_database: HashMap<u64, Vec<u64>> = HashMap::new();
    let all_operators: Vec<Vec<(usize, usize)>> = topology.get_all_operators();
    let mut initial = MatrixF2::identity_half(&topology.nvertices, &topology.nvertices);
    let mut depth = 0;
    global_database.insert(initial.cannonical_representation_step1(), (0, 0));
    target_database.insert(initial.cannonical_representation_step1(), Vec::new());
    let mut queue: Vec<(MatrixF2, u64)> = vec![(initial, 0)];
    loop {
        depth += 1;
        if queue.is_empty() {
            break;
        }
        let mut new_queue: Vec<(MatrixF2, u64)> = Vec::new();
        for (mut matrix, _) in queue.drain(..) {
            let new_parent = matrix.cannonical_representation_step1();
            for operator in all_operators.iter() {
                matrix.apply_operator(operator);
                let repr = matrix.cannonical_representation_step1();
                if !global_database.contains_key(&repr) {
                    global_database.insert(repr, (depth, new_parent));
                    new_queue.push((matrix.clone(), new_parent));
                    let genealogy = reconstruct_genealogy(&global_database, &repr);
                    target_database.insert(repr, genealogy);
                }
                matrix.apply_operator(operator);
            }
        }
        queue = new_queue;
    }
    return target_database;
}

pub fn enumerate_step_3(topology: &Topology) -> HashMap<u64, Vec<u64>> {
    let topology = topology.get_block_topology();
    let all_operators = topology.get_all_operators();
    let mut global_database: HashMap<u64, (i32, u64)> = HashMap::new();
    let mut target_database: HashMap<u64, Vec<u64>> = HashMap::new();
    let initial = MatrixF2::identity(&topology.nvertices, &topology.nvertices);
    let mut depth = 0;
    global_database.insert(initial.cannonical_representation_step3(), (0, 0));
    target_database.insert(initial.cannonical_representation_step3(), Vec::new());
    let mut queue: Vec<(MatrixF2, u64)> = vec![(initial, 0)];
    loop {
        depth += 1;
        if queue.is_empty() {
            break;
        }
        let mut new_queue: Vec<(MatrixF2, u64)> = Vec::new();
        for (mut matrix, _) in queue.drain(..) {
            let new_parent = matrix.cannonical_representation_step3();
            for operator in all_operators.iter() {
                matrix.apply_operator(operator);
                let repr = matrix.cannonical_representation_step3();
                if !global_database.contains_key(&repr) {
                    global_database.insert(repr, (depth, new_parent));
                    new_queue.push((matrix.clone(), new_parent));
                    let genealogy = reconstruct_genealogy(&global_database, &repr);
                    target_database.insert(repr, genealogy);
                }
                matrix.apply_operator(operator);
            }
        }
        queue = new_queue;
    }
    return target_database;
}

/// Generates the step2 database for a given topology and dumps the resulting database in some file
pub fn generate_full_database_step_2(
    topology: &Topology,
    output_fname: &str,
) -> std::io::Result<()> {
    let database = enumerate_step_2(topology);
    let all_operators = topology.get_all_operators();
    let mut file = File::create(output_fname)?;
    for repr in database.keys() {
        writeln!(
            file,
            "{} {:?}",
            repr,
            get_circuit_from_repr(
                &database,
                &all_operators,
                repr,
                topology.nvertices,
                topology.nvertices
            )
        )?;
    }
    return Result::Ok(());
}

pub fn generate_full_database_step_3(
    topology: &Topology,
    output_fname: &str,
) -> std::io::Result<()> {
    let database = enumerate_step_3(topology);
    let all_operators = topology.get_block_topology().get_all_operators();
    let mut file = File::create(output_fname)?;
    for repr in database.keys() {
        writeln!(
            file,
            "{} {:?}",
            repr,
            get_circuit_from_repr(
                &database,
                &all_operators,
                repr,
                topology.nvertices / 2,
                topology.nvertices / 2
            )
        )?;
    }
    return Result::Ok(());
}

/// Generates any database for a given topology and dumps the resulting database in some file
#[pyfunction]
pub fn generate_database(
    topology_name: &str,
    output_fname: &str,
    step: i32,
) -> std::io::Result<()> {
    let topology = Topology::from_string(topology_name);
    let database = match step {
        1 => enumerate_step_1(&topology),
        2 => enumerate_step_2(&topology),
        3 => enumerate_step_3(&topology),
        _ => panic!("There is not step {}", step),
    };
    let all_operators = if step < 3 {
        topology.get_all_operators()
    } else {
        topology.get_block_topology().get_all_operators()
    };
    let mut file = File::create(output_fname)?;
    for repr in database.keys() {
        writeln!(
            file,
            "{},{}",
            repr,
            get_circuit_from_repr(
                &database,
                &all_operators,
                repr,
                if step < 3 {
                    topology.nvertices
                } else {
                    topology.nvertices / 2
                },
                if step < 3 {
                    topology.nvertices
                } else {
                    topology.nvertices / 2
                },
            )
            .iter()
            .map(|(a, b)| format!("{a}-{b}"))
            .collect::<Vec<String>>()
            .join("|")
        )?;
    }
    return Result::Ok(());
}

/// Loads a database out of a file name
#[pyfunction]
pub fn load_database(fname: &str) -> HashMap<u64, Vec<(usize, usize)>> {
    let file = File::open(fname).unwrap();
    let lines = BufReader::new(file).lines();
    let mut database = HashMap::new();
    for line in lines {
        if let Result::Ok(line) = line {
            let mut line = line.split(',');
            let repr = line.next().unwrap().parse::<u64>().unwrap();
            if let Some(cnots) = line.next() {
                if cnots.len() > 1 {
                    let cnots: Vec<(usize, usize)> = cnots
                        .split('|')
                        .map(|s| {
                            let mut blob = s.split('-');
                            let a = blob.next().unwrap().parse::<usize>().unwrap();
                            let b = blob.next().unwrap().parse::<usize>().unwrap();
                            (a, b)
                        })
                        .collect();
                    database.insert(repr, cnots);
                } else {
                    database.insert(repr, Vec::new());
                }
            } else {
                database.insert(repr, Vec::new());
            }
        }
    }
    return database;
}

#[pyfunction]
/// The matrix is assumed to be in col major
pub fn get_circuit_from_matrix(
    matrix: Vec<u8>,
    database: HashMap<u64, Vec<(usize, usize)>>,
    step: i32,
) -> Vec<(usize, usize)> {
    let n = if step > 1 {
        (matrix.len() as f32).sqrt() as usize
    } else {
        2 * ((matrix.len() / 2) as f64).sqrt() as usize
    };
    let m = if step > 1 { n } else { n / 2 };
    let mut mat = MatrixF2::zero(&n, &n);
    for i in 0..n {
        for j in 0..m {
            if matrix[j * n + i] == 1 {
                mat.set(i, j);
            }
        }
    }
    return database[&mat.cannonical_representation()].clone();
}

#[pymodule]
fn block_kutin(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(enumeration::generate_database))?;
    m.add_wrapped(wrap_pyfunction!(enumeration::load_database))?;
    m.add_wrapped(wrap_pyfunction!(enumeration::get_circuit_from_matrix))?;
    Ok(())
}
