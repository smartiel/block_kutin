use super::matrix::MatrixF2;
use super::operators::Topology;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;

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
    println!(
        "There are {} depth-1 operators for this topology ({} matchings)",
        all_operators.len(),
        topology.get_all_matchings().len()
    );
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
        println!("Enumerating for depth {depth}");
        println!("There are {} operators of depth {}", queue.len(), depth - 1);
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
    println!(
        "There are {} depth-1 operators for this topology ({} matchings)",
        all_operators.len(),
        topology.get_all_matchings().len()
    );
    let m = topology.nvertices / 2;
    let mut initial = MatrixF2::identity(&topology.nvertices, &m);
    let mut depth = 0;
    global_database.insert(initial.cannonical_representation_step1(), (0, 0));
    let mut queue: Vec<(MatrixF2, u64)> = vec![(initial, 0)];
    loop {
        depth += 1;
        println!("Enumerating for depth {depth}");
        println!("There are {} operators of depth {}", queue.len(), depth - 1);
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

/// Generates the step1 database for a given topology and dumps the resulting database in some file
pub fn generate_full_database_step_1(
    topology: &Topology,
    output_fname: &str,
) -> std::io::Result<()> {
    let database = enumerate_step_1(topology);
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
                topology.nvertices / 2
            )
        )?;
    }
    return Result::Ok(());
}
